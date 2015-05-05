#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use CortexCovgFile;

#
# Read in and print out .colour_covg cortex output
#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 13 Nov 2011
#

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./test_read_align.pl [.colour_covgs]
  Example script for reading in alignments from .colour_covg files\n";

  exit;
}

if(@ARGV > 1)
{
  print_usage();
}

my $covg_file = shift;

#
# Open .colour_covgs handle
#
my $covg_handle;

if(defined($covg_file) && $covg_file ne "-")
{
  open($covg_handle, $covg_file)
    or print_usage("Cannot open .colour_covgs file '$covg_file'");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covg_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a Cortex .colour_covgs file");
}

#
# Read .colour_covg data
#

my $covgfile = new CortexCovgFile($covg_handle);

# Start reading aligned entries
my ($read_name, $sequence, $colours_arrref);

while((($read_name, $sequence, $colours_arrref) = $covgfile->read_align_entry()) &&
      defined($read_name))
{
  print_colour_covg_entry($read_name, $sequence, $colours_arrref);
}

close($covg_handle);

