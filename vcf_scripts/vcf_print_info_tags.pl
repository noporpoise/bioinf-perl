#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;

# Config #
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_info_tags.pl <file.vcf> <infotag1 ..>\n";
  print STDERR "  Prints comma separated info tag values from VCF entries\n";
  print STDERR "  If <file.vcf> is '-' reads from STDIN\n";
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

my $vcf_file = shift;

my @tags = @ARGV;

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

# Print CSV header
print join($csvsep, @tags)."\n";

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my @data = ();

  for my $tag (@tags)
  {
    my $d = defined($vcf_entry->{'INFO'}->{$tag}) ? $vcf_entry->{'INFO'}->{$tag} : "";
    push(@data, $d);
  }

  print join($csvsep, @data)."\n";
}

$vcf->vcf_close();
