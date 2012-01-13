#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_header.pl [--entries] [vcf] [extra_header]\n";
  print STDERR "  Print vcf file header (reads from STDIN if <file> is '-').\n";
  print STDERR "  Adds extra headers if passed after file.\n";
  print STDERR "\n";
  print STDERR "  [extra_header] should be of the following form:\n";
  print STDERR "    INFO,ID,Number,Type,Description    e.g. 'INFO,T,1,String,T number'\n";
  print STDERR "    FORMAT,ID,Number,Type,Description  e.g. 'FORMAT,T,1,String,T number'\n";
  print STDERR "    FILTER,ID,Description              e.g. 'FILTER,T,T number'\n";
  print STDERR "    ALT,ID,Description                 e.g. 'ALT,T,T number'\n";
  print STDERR "\n";
  print STDERR "  --entries  Print VCF entries as well as header\n";
  exit;
}

my $print_entries = 0;

if(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-e(ntries)?$/i)
  {
    $print_entries = 1;
    shift;
  }
  elsif($ARGV[0] =~ /^-?-h(elp)?$/i)
  {
    print_usage();
  }
}

my $vcf_file = shift;

my @extra_headers = @ARGV;

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

for my $extra_header (@extra_headers)
{
  my @parts = split(",", $extra_header);

  my ($tag_column, $tag_id, $tag_number, $tag_type, $tag_description);

  if(@parts == 3)
  {
    ($tag_column,$tag_id,$tag_description) = @parts;
  }
  elsif(@parts == 5)
  {
    ($tag_column, $tag_id, $tag_number, $tag_type, $tag_description) = @parts;
  }
  else
  {
    print_usage("Invalid extra header argument: '$extra_header'");
  }

  if($tag_description =~ /^\".*\"$/)
  {
    # Trim first and last characters off
    $tag_description = substr($tag_description, 1, -1);
  }
  elsif($tag_description =~ /\"/)
  {
    print_usage("Description should not contain quotation marks");
  }

  $vcf->add_header_tag($tag_column, $tag_id, $tag_number,
                       $tag_type, $tag_description);
}

$vcf->print_header();

if($print_entries)
{
  my $vcf_entry;

  while(defined($vcf_entry = $vcf->read_entry()))
  {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
