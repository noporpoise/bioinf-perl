#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use UsefulModule; # for num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_filter_indels.pl [--clean|--invert] [in.vcf]\n";
  print STDERR "  Prints variants were SVLEN != 0\n";
  print STDERR "  --clean           One allele is 0bp and the other is >0bp\n";
  print STDERR "  --invert          Print non-indels\n";
  print STDERR "  --clean --invert  Print non-clean indels\n";
  exit;
}

if(@ARGV > 3)
{
  print_usage();
}

my $filter_clean = 0;
my $filter_invert = 0;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-c(lean)?$/i)
  {
    $filter_clean = 1;
    shift;
  }
  elsif($ARGV[0] =~ /^-?-inv(ert)?$/i)
  {
    $filter_invert = 1;
    shift;
  }
  else
  {
    last;
  }
}

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
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

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $is_indel = $vcf_entry->{'INFO'}->{'SVLEN'} != 0;
  my $one_is_zero = (length($vcf_entry->{'true_REF'}) == 0 ||
                     length($vcf_entry->{'true_ALT'}) == 0);

  my $print = $is_indel && (!$filter_clean || $one_is_zero);

  if($print != $filter_invert)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_printed++;
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_of_printed / $num_of_variants;

print STDERR "vcf_filter_indels.pl: " .
             num2str($num_of_printed) . " / " . num2str($num_of_variants) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
