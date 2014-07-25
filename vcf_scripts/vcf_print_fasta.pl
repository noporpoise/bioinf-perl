#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

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
    --pad <k>        Pad alleles with k bases of surrounding sequence
    --trim <t>       Trim sequences longer than <t>\n";

  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 2)
{
  print_usage();
}

my $trim;
my $pad_alleles = 0;
my $use_ancestral = 0;
my $split_fasta = 0;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-trim$/i)
  {
    shift;
    $trim = shift;

    if($trim !~ /^\d+$/ || $trim == 0)
    {
      print_usage("Invalid trim value ('$trim') - " .
                  "must be positive integer (>0)");
    }
  }
  elsif($ARGV[0] =~ /^-?-pad$/i)
  {
    shift;
    $pad_alleles = shift;

    if($pad_alleles !~ /^\d+$/ || $pad_alleles == 0)
    {
      print_usage("Invalid --pad value ('$pad_alleles') - " .
                  "must be positive integer (>=0)");
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

if(@ref_files == 0)
{
  push(@ref_files, '-');
}

if(!defined($vcf_file))
{
  print_usage();
}

if($max_flank_size !~ /^\d+$/)
{
  print_usage("Max flank size must be a positive integer (>=0): $max_flank_size");
}
elsif(defined($trim) && $trim > $max_flank_size)
{
  print_usage("--trim cannot be larger than max flank size " .
              "(trim: $trim; flank: $max_flank_size)");
}

$max_flank_size += $pad_alleles;

#
# Read VCF
#
my $vcf = vcf_open($vcf_file);

#
# Load reference
#
my $genome = new RefGenome();
$genome->load_from_files(@ref_files);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Get flanks
  my $lflank = "";
  my $rflank = "";

  if($max_flank_size > 0)
  {
    if(@ref_files == 0)
    {
      if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
         !defined($vcf_entry->{'INFO'}->{'right_flank'}))
      {
        print_usage("Flanks missing on entry $vcf_entry->{'ID'}");
      }

      $lflank = $vcf_entry->{'INFO'}->{'left_flank'};
      $rflank = $vcf_entry->{'INFO'}->{'right_flank'};

      if(length($lflank) > $max_flank_size)
      {
        $lflank = substr($lflank, -$max_flank_size);
      }

      if(length($rflank) > $max_flank_size)
      {
        $rflank = substr($rflank, 0, $max_flank_size);
      }
    }
    else
    {
      ($lflank, $rflank) = vcf_get_flanks($vcf_entry, $genome, $max_flank_size);
    }
  }

  # Get alleles
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};

  if($pad_alleles > 0)
  {
    # Work out how much padding we can get either side
    my $lpad = min(length($lflank), $pad_alleles);
    my $rpad = min(length($rflank), $pad_alleles);
    # Add padding to the alleles
    $ref_allele = substr($lflank, -$lpad).$ref_allele.substr($rflank, 0, $rpad);
    $alt_allele = substr($lflank, -$lpad).$alt_allele.substr($rflank, 0, $rpad);
    # Remove bases used in padding from the flanks
    $lflank = substr($lflank, 0, length($lflank)-$lpad);
    $rflank = substr($rflank, $rpad);
  }

  # Get alleles
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
    print_to_fasta($vcf_entry->{'ID'}."_flank5", $lflank);
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_anc" : "_ref"), $allele1);
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_der" : "_alt"), $allele2);
    print_to_fasta($vcf_entry->{'ID'}."_flank3", $rflank);
  }
  else
  {
    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_anc" : "_ref"),
                   $lflank . $allele1 . $rflank);

    print_to_fasta($vcf_entry->{'ID'}.($using_ancestral ? "_der" : "_alt"),
                   $lflank . $allele2 . $rflank);
  }
}

$vcf->vcf_close();

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
