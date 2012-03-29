#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

## Config
my $align_cmd = "needleman_wunsch --pretty ";
##

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_view.pl [OPTIONS] [in.vcf]

OPTIONS:
  --snps      Variants must have snps
  --indels    Variants must have indels
  --no_snps   Variants must not have snps
  --no_indels Variants must not have indels
  --as_vcf    Print in VCF format\n";

  exit;
}

if(@ARGV > 4)
{
  print_usage();
}

my $require_snps = 0;
my $no_snps = 0;

my $require_indels = 0;
my $no_indels = 0;

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
  else
  {
    last;
  }
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
  }
  elsif(defined(get_clean_indel($vcf_entry)))
  {
    $has_indels = 1;

    if(length($allele1) == 0)
    {
      $allele1 = "-" x length($allele2);
      $sep = " " x length($allele2);
    }
    else
    {
      $allele2 = "-" x length($allele1);
      $sep = " " x length($allele1);
    }
  }
  else
  {
    my $alignment = `$align_cmd $allele1 $allele2`;
    ($allele1, $sep, $allele2) = split("\n", $alignment);

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
      print $vcf_entry->{'ID'}.":\n";
      print $allele1."\n";
      #print $sep."\n";
      print $allele2."\n";
      print "\n";
    }
  }
}

close($vcf_handle);
