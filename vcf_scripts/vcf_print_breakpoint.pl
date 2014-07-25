#!/usr/bin/perl

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
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_breakpoint.pl <flank_size> [in.vcf]\n";

  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $flank_size = shift;
my $vcf_file = shift;

if($flank_size !~ /^\d+$/ || $flank_size <= 0)
{
  print_usage("Invalid flank size '$flank_size'");
}

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  if(defined($aa) && ($aa eq "0" || $aa eq "1"))
  {
    my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
    my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

    if(!defined($left_flank) || !defined($right_flank))
    {
      print_usage("Missing left/right flank '".$vcf_entry->{'ID'}."'");
    }
    elsif(length($left_flank) < $flank_size || length($right_flank) < $flank_size)
    {
      print_usage("Flank is too short on variant '".$vcf_entry->{'ID'}."'");
    }

    print ">".$vcf_entry->{'ID'}."\n";
    print substr($left_flank, -$flank_size) .
          $vcf_entry->{$aa eq "0" ? 'true_REF': 'true_ALT'} .
          substr($right_flank, 0, $flank_size) . "\n";
  }
}

$vcf->vcf_close();
