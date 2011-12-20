#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_header.pl [--entries] [vcf] [extra_header1 ...]\n";
  print STDERR "  Print vcf file header (reads from STDIN if <file> is '-').\n";
  print STDERR "  Adds extra headers if passed after file.\n";
  print STDERR "  -entries => print VCF entries as well as header\n";
  exit;
}

my $print_entries = 0;

if(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-e(ntries)?/i)
  {
    $print_entries = 1;
    shift;
  }
  elsif($ARGV[0] =~ /-?-h(elp)?/i)
  {
    print_usage();
  }
}

my $vcf_file = shift;

my @extra_header_lines = @ARGV;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

my $vcf = new VCFFile($vcf_handle);

print vcf_add_to_header($vcf->get_header(), @extra_header_lines);

if($print_entries)
{
  my $vcf_entry;

  while(defined($vcf_entry = $vcf->read_entry()))
  {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
