#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use CortexCovgFile;

#
# Filter novel sequences based on covg
#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 25 Jun 2013
#

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "" .
"Usage: ./cortex_filter_novel.pl [--filter <colour> <covg>] [.colour_covgs]
  Filter novel sequences based on covg\n";
  exit(-1);
}

my $covg_file;
my @filter_cols = ();
my @filter_covgs = ();

while(@ARGV > 0)
{
  my $arg = shift;
  if(@ARGV == 0) { $covg_file = $arg; last; }
  elsif($arg =~ /^--filter$/i) {
    my $col = shift;
    my $covg = shift;
    if(!defined($col) || $col !~ /^\d+$/ || !defined($covg) || $covg !~ /^\d+$/) {
      print_usage("Invalid --filter arguments");
    }
    push(@filter_cols, $col);
    push(@filter_covgs, $covg);
  }
  else { print_usage("Unknown option: $arg"); }
}

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
  if(get_filter_status($colours_arrref)) {
    print_colour_covg_entry($read_name, $sequence, $colours_arrref);
  }
}

close($covg_handle);

sub get_filter_status
{
  my ($covgs) = @_;
  my $len = @{$covgs->[0]};
  my $num_filters = @filter_cols;
  if($num_filters == 0 || $len <= 2) { return 1; }

  for(my $i = 0; $i < $num_filters; $i++) {
    my ($col, $covg) = ($filter_cols[$i], $filter_covgs[$i]);
    if($col >= @$covgs) { die("colour $col > ".(@$covgs-1)."\n"); }
    if(get_covg_median($covgs->[$col]) >= $covg) { return 1; }
  }
  return 0;
}

sub get_covg_median
{
  my ($arr) = @_;
  my $len = @$arr;
  my @tmp = sort map {$arr->[$_]} (1..($len-2));
  if(@tmp == 0) { return 0; }
  if(@tmp % 2 == 0) { return ($tmp[$#tmp/2] + $tmp[$#tmp/2+1]) / 2; }
  else { return $tmp[$#tmp/2]; }
}
