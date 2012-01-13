#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use UsefulModule; # num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./vcf_filter_by_checks.pl [--pass|--fail] [in.vcf]\n";
  print STDERR "  Prints variants that passed (or failed) the filtering\n";
  print STDERR "  If [in.vcf] is '-', reads from stdin\n";
  exit;
}

if(@ARGV > 2) {
  print_usage();
}

my $print_passed = 1;
my $vcf_file;

if(@ARGV == 2) {
  my $arg_cmd = shift;
  $vcf_file = shift;

  if($arg_cmd =~ /^-?-p(ass)?$/i) {
    $print_passed = 1;
  }
  elsif($arg_cmd =~ /^-?-f(ail)?$/i) {
    $print_passed = 0;
  }
  else {
    print_usage();
  }
}
else {
  $vcf_file = shift;
}

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

my $num_of_filtered_entries = 0;
my $total_num_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry))
{
  $total_num_entries++;

  if(($vcf_entry->{'FILTER'} =~ /^PASS$/i) == $print_passed)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_of_filtered_entries / $total_num_entries;

print STDERR "vcf_filter_by_checks.pl: " . num2str($num_of_filtered_entries) .
             " / " . num2str($total_num_entries) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
