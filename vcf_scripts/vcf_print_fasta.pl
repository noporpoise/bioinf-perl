#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use RefGenome;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_print_fasta.pl [OPTIONS] <flank_size> <in.vcf> [ref.files ..]
  Prints alleles plus flanking sequence from VCF file in FASTA format. If
  <in.vcf> is '-', reads from stdin.  If [ref.files] omitted, uses INFO flank
  tags (left_flank=.. & right_flank=..).

  OTPIONS:
    --ancestral      Use AA tags to print ancestral,derived instead of ref,alt
                     When AA tags not available, reverts to ref,alt
    --split          Prints 5' flank, alleles and 3' flank on different lines
    --trim <t>       Trim sequences longer than <t>\n";

  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 2 || @ARGV > 5)
{
  print_usage();
}

my $trim;
my $use_ancestral = 0;
my $split_fasta = 0;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-trim$/i)
  {
    $trim = shift(@ARGV);

    if($trim !~ /^\d+$/ || $trim == 0)
    {
      print_usage("Invalid trim value ('$trim') - " .
                  "must be positive integer (>0)");
    }
  }
  elsif($ARGV[0] =~ /^-?-ancestral$/i)
  {
    shift;
    $use_ancestral = 1;
  }
  elsif($ARGV[0] =~ /^-?-split$/i)
  {
    shift;
    $split_fasta = 1;
  }
  elsif($ARGV[0] =~ /^-/)
  {
    print_usage("Unknown argument '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

my $max_flank_size = shift;
my $vcf_file = shift;
my @ref_files = @ARGV;

if($max_flank_size !~ /^\d+$/)
{
  print_usage("Max flank size must be a positive integer (>=0): $max_flank_size");
}
elsif(defined($trim) && $trim > $max_flank_size)
{
  print_usage("--trim cannot be larger than max flank size " .
              "(trim: $trim; flank: $max_flank_size)");
}

my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

my $genome = new RefGenome();

if(@ref_files > 0)
{
  # Load reference
  $genome->load_from_files(@ref_files);
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Get flanks
  my $left_flank = "";
  my $right_flank = "";

  if($max_flank_size > 0)
  {
    if(@ref_files == 0)
    {
      if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
         !defined($vcf_entry->{'INFO'}->{'right_flank'}))
      {
        print_usage("Flanks missing on entry $vcf_entry->{'ID'}");
      }

      $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
      $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};
  
      if(length($left_flank) > $max_flank_size)
      {
        $left_flank = substr($left_flank, 0, $max_flank_size);
      }
    
      if(length($left_flank) > $max_flank_size)
      {
        $right_flank = substr($right_flank, -$max_flank_size);
      }
    }
    else
    {
      ($left_flank, $right_flank) = vcf_get_flanks($vcf_entry, $genome, $max_flank_size);
    }
  }

  # Get alleles
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};
  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  my ($allele1, $allele2) = ($ref_allele, $alt_allele);

  my $using_ancestral = 0;
  if($use_ancestral && defined($aa) && ($aa == 0 || $aa == 1))
  {
    $using_ancestral = 1;

    if($aa == 1)
    {
      ($allele1, $allele2) = ($alt_allele, $ref_allele);
    }# else already correct
  }

  if($split_fasta)
  {
    print_to_fasta($vcf_entry->{'ID'}."_flank5", $left_flank);
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_anc" : "_ref"), $allele1);
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_der" : "_alt"), $allele2);
    print_to_fasta($vcf_entry->{'ID'}."_flank3", $right_flank);
  }
  else
  {
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_anc" : "_ref"),
                   $left_flank . $allele1 . $right_flank);

    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_der" : "_alt"),
                   $left_flank . $allele2 . $right_flank);
  }
}

close($vcf_handle);

sub print_to_fasta
{
  my ($fasta_name,$fasta_seq) = @_;

  print ">$fasta_name\n";
  
  if(defined($trim))
  {
    print substr($fasta_seq, 0, $trim) . "\n";
  }
  else
  {
    print "$fasta_seq\n";
  }
}
