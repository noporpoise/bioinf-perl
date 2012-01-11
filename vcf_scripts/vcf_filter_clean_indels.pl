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

  print STDERR "Usage: ./vcf_filter_clean_indels.pl [--invert] [in.vcf]\n";
  print STDERR "  Prints variants were one allele is 0bp and the other is >0bp\n";
  exit;
}

if(@ARGV > 2)
{
  print_usage();
}

my $invert = 0;

if(@ARGV > 0 && $ARGV[0] =~ /^-?-invert?/i)
{
  $invert = 1;
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

print $vcf->get_header();

my $num_of_variants = 0;
my $num_printed = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $is_clean = ($vcf_entry->{'INFO'}->{'SVLEN'} > 0 &&
                  (length($vcf_entry->{'true_REF'}) == 0 ||
                   length($vcf_entry->{'true_ALT'}) == 0));

  if($is_clean != $invert)
  {
    $num_printed++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_printed / $num_of_variants;

print STDERR "vcf_get_clean_indels.pl: " .
             num2str($num_printed) . " / " . num2str($num_of_variants) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
