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
"Usage: ./vcf_view.pl [--only_snps|--only_indels|--only_both|--as_vcf] [in.vcf]
  --only_snps   Only print SNPs
  --only_indels Only print 
  --only_both   \n";

  exit;
}

if(@ARGV > 3)
{
  print_usage();
}

my $only_snps = 0;
my $only_indels = 0;
my $only_both = 0;
my $view_all = 1;

my $view_as_vcf = 0;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-only_snps$/i)
  {
    shift;
    $only_snps = 1;
    $view_all = 0;
  }
  elsif($ARGV[0] =~ /^-?-only_indels$/i)
  {
    shift;
    $only_indels = 1;
    $view_all = 0;
  }
  elsif($ARGV[0] =~ /^-?-only_both$/i)
  {
    shift;
    $only_both = 1;
    $view_all = 0;
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

if(sum($only_snps, $only_indels, $only_both, $view_all) != 1)
{
  print_usage("Can only specify one of --only_snps, --only_indels, --only_both");
}

my $vcf_file = shift;

if(@ARGV > 0)
{
  print_usage("Unexpected option '$ARGV[0]'");
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
  if($only_snps) {
    $vcf->add_header_metainfo("vcf_view", "Only SNPs");
  }
  elsif($only_indels) {
    $vcf->add_header_metainfo("vcf_view", "Only indels");
  }
  elsif($only_both) {
    $vcf->add_header_metainfo("vcf_view", "Only both");
  }

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

  if($view_all ||
     $only_snps && $has_snps && !$has_indels ||
     $only_indels && !$has_snps && $has_indels ||
     $only_both && $has_snps && $has_indels)
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
