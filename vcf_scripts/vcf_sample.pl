#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_sample.pl <prob> [file.vcf]
  Sample from a VCF file, taking each variant with a probability of <prob>.  
  must be: 0 <= Prob <= 1\n";

  exit(-1);
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $prob = shift;
my $vcf_file = shift;

if($prob !~ /^(\d*\.\d+|\d+\.?\d*)$/ || $prob > 1 || $prob < 0)
{
  print_usage("Invalid probability '$prob' (must be: 0 <= Prob <= 1)");
}

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  if(rand() < $prob)
  {
    $vcf->print_entry($vcf_entry);
  }
}

$vcf->vcf_close();
