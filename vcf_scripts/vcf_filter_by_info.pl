#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_filter_by_info.pl <file.vcf> " .
               "[<INFO_FIELD> <VALUE_REGEX>] ..\n";
  print STDERR "  Isaac Turner <isaac.turner\@dtc.ox.ac.uk> 2011/03/26\n";
  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $vcf_file = shift;

if((@ARGV % 2) != 0)
{
  print_usage();
}

my %searches = ();

for(my $i = 0; $i < @ARGV; $i+=2)
{
  $searches{$ARGV[$i]} = $ARGV[$i+1];
}

#
# Open VCF Handle
#
my $vcf_handle;

if($vcf_file ne "-") {
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

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $info_hashref = $vcf_entry->{'INFO'};
  my $flags_hashref = $vcf_entry->{'INFO_flags'};

  my $match = 1;

  my ($key,$search);

  for my $key (keys %searches)
  {
    if(!defined($flags_hashref->{$key}) &&
       (!defined($info_hashref->{$key}) ||
        $info_hashref->{$key} !~ /$searches{$key}/i))
    {
      $match = 0;
      last;
    }
  }

  if($match == 1) {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
