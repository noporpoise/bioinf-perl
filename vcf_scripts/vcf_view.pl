#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use FASTNFile;

## Config
my $align_cmd = "needleman_wunsch --pretty";
my $mismatch_colour = "\033[92m"; # Mismatch (GREEN)
my $indel_colour = "\033[91m"; # Insertion / deletion (RED)
my $stop_colour = "\033[0m";
##

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_view.pl [OPTIONS] [in.vcf]

OPTIONS:
  --snps       Variants must have snps
  --indels     Variants must have indels
  --no_snps    Variants must not have snps
  --no_indels  Variants must not have indels
  
  
  --ref <ref>  Load reference genome
  --flanks <f> Print <f> bp of flank, using left_flank INFO tag or --ref file
  --colour     Print with colour
  
  --as_vcf     Print in VCF format\n";

  exit;
}

my $require_snps = 0;
my $no_snps = 0;

my $require_indels = 0;
my $no_indels = 0;

my $ref_file;
my $flank_size = 0;

my $colour_output = 0;
my $view_as_vcf = 0;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-snps$/i)
  {
    shift;
    $require_snps = 1;
  }
  elsif($ARGV[0] =~ /^-?-no_snps$/i)
  {
    shift;
    $no_snps = 1;
  }
  elsif($ARGV[0] =~ /^-?-indels$/i)
  {
    shift;
    $require_indels = 1;
  }
  elsif($ARGV[0] =~ /^-?-no_indels$/i)
  {
    shift;
    $no_indels = 1;
  }
  elsif($ARGV[0] =~ /^-?-as_vcf$/i)
  {
    shift;
    $view_as_vcf = 1;
  }
  elsif($ARGV[0] =~ /^-?-ref$/i)
  {
    shift;
    $ref_file = shift;
    
    if(!defined($ref_file))
    {
      print_usage("missing --ref <ref> argument");
    }
  }
  elsif($ARGV[0] =~ /^-?-flanks?$/i)
  {
    shift;
    $flank_size = shift;
    
    if(!defined($flank_size) || $flank_size !~ /^\d+$/)
    {
      print_usage("--flanks requires a positive integer argument");
    }
  }
  elsif($ARGV[0] =~ /^-?-colour$/i)
  {
    shift;
    $colour_output = 1;
    $align_cmd .= " --colour";
  }
  else
  {
    last;
  }
}

if($colour_output && $view_as_vcf)
{
  print_usage("Cannot use --as_vcf and --colour together");
}

if($require_snps && $no_snps)
{
  print_usage("Cannot specify both --no_snps and --snps");
}

if($require_indels && $no_indels)
{
  print_usage("Cannot specify both --no_indels and --indels");
}

my $vcf_file = shift;

if(@ARGV > 0)
{
  print_usage("excess option '$ARGV[0]'");
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Load reference
#
my %ref_genomes;

if(defined($ref_file))
{
  my ($ref_genomes_hashref) = read_all_from_files($ref_file);
  %ref_genomes = %$ref_genomes_hashref;
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

if($view_as_vcf)
{
  my $description = "";

  if($no_snps)
  {
    $description .= "NoSNPs";
  }
  elsif($require_snps)
  {
    $description .= "RequireSNPs";
  }
  
  if($no_snps)
  {
    $description .= "NoIndels";
  }
  elsif($require_snps)
  {
    $description .= "RequireIndels";
  }
  
  if($description eq "")
  {
    $description = "All";
  }

  $vcf->add_header_metainfo("vcf_view", $description);
  $vcf->print_header();
}

my $vcf_entry;
my ($allele1, $sep, $allele2);

while(defined($vcf_entry = $vcf->read_entry()))
{
  $allele1 = $vcf_entry->{'true_REF'};
  $allele2 = $vcf_entry->{'true_ALT'};

  if(defined($vcf_entry->{'INFO'}->{'AA'}) &&
     $vcf_entry->{'INFO'}->{'AA'} eq "1")
  {
    ($allele1, $allele2) = ($allele2, $allele1);
  }

  my $has_snps = 0;
  my $has_indels = 0;

  if(is_snp($vcf_entry))
  {
    $has_snps = 1;
    $sep = "*";
    
    if($colour_output)
    {
      $allele1 = $mismatch_colour.$allele1.$stop_colour;
      $allele2 = $mismatch_colour.$allele2.$stop_colour;
    }
  }
  elsif(defined(get_clean_indel($vcf_entry)))
  {
    $has_indels = 1;

    if(length($allele1) == 0)
    {
      $allele1 = "-" x length($allele2);
      $sep = " " x length($allele2);
      
      if($colour_output) {
        $allele2 = $indel_colour.$allele2.$stop_colour;
      }
    }
    else
    {
      $allele2 = "-" x length($allele1);
      $sep = " " x length($allele1);
      
      if($colour_output) {
        $allele1 = $indel_colour.$allele1.$stop_colour;
      }
    }
  }
  else
  {
    my $alignment = `$align_cmd $allele1 $allele2`;
    ($allele1, $sep, $allele2) = split("\n", $alignment);

    # To remove colouring
    #my $cpy1 = $allele1;
    #$cpy1 =~ s/\\033[\d+m//gi;

    $has_snps = ($sep =~ /\*/);
    $has_indels = ($sep =~ / /);
  }

  if(!($no_snps && $has_snps) &&
     !($no_indels && $has_indels) &&
     (!$require_snps || $has_snps) &&
     (!$require_indels || $has_indels))
  {
    if($view_as_vcf)
    {
      $vcf->print_entry($vcf_entry);
    }
    else
    {
      my $left_flank = "";
      my $right_flank = "";

      if($flank_size > 0)
      {
        my $chr = $vcf_entry->{'CHROM'};

        $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
        $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

        if(defined($left_flank) && length($left_flank) >= $flank_size &&
           defined($right_flank) && length($right_flank) >= $flank_size)
        {
          $left_flank = substr($left_flank, -$flank_size);
          $right_flank = substr($right_flank, 0, $flank_size);

          $sep = ("|" x $flank_size) . $sep . ("|" x $flank_size);
        }
        elsif(defined($ref_file) && defined($ref_genomes{$chr}))
        {
          $left_flank = substr($ref_genomes{$chr},
                               $vcf_entry->{'true_POS'} - 1 - $flank_size,
                               $flank_size);

          $right_flank = substr($ref_genomes{$chr},
                                $vcf_entry->{'true_POS'} - 1 +
                                  length($vcf_entry->{'true_REF'}),
                                $flank_size);

          $sep = ("|"x$flank_size) . $sep . ("|"x$flank_size);
        }
      }

      print $vcf_entry->{'ID'}.":\n";
      print $left_flank.$allele1.$right_flank."\n";
      print $sep."\n";
      print $left_flank.$allele2.$right_flank."\n";
      print "\n";
    }
  }
}

close($vcf_handle);
