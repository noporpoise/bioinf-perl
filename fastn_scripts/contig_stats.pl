#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "" .
"Usage: ./contig_stats.pl [--print-csv] <in.fa|fq>
  Length distribution stats.
";

  exit(-1);
}

my $print_csv = 0;

while(@ARGV > 1) {
  if($ARGV[0] =~ /^--print-csv$/) {
    shift(@ARGV);
    $print_csv = 1;
  } else {
    print_usage("Bad argument: $ARGV[0]");
  }
}

if(@ARGV != 1) { print_usage(); }

my $path = shift(@ARGV);
my @lengths = ();

my $fastn = open_fastn_file($path);
my ($title,$seq);
while((($title,$seq) = $fastn->read_next()) && defined($title)) {
  push(@lengths, length($seq));
}
close_fastn_file($fastn);

# min max median mode n50
@lengths = sort{$a <=> $b} @lengths;
my $ncontigs = scalar(@lengths);
my $sum = sum(@lengths);

if($ncontigs == 0) { print STDERR "[contig_stats.pl] No sequences\n"; exit -1; }

my $median = find_median(@lengths);
my $mode = find_mode(@lengths);
my $n50 = find_N50($sum,@lengths);

# Printing
my $lwidth = 7;
my $rwidth = 15;

if($print_csv) { print "metric,value\n"; }

# Some lines $linewidth+2 for 1 decimal place
print_cols("contigs", $lwidth, num2str($ncontigs),              $rwidth);
print_cols("length",  $lwidth, num2str($sum),                   $rwidth);
print_cols("min",     $lwidth, num2str($lengths[0]),            $rwidth);
print_cols("max",     $lwidth, num2str($lengths[$ncontigs-1]),  $rwidth);
print_cols("mean",    $lwidth, num2str($sum/$ncontigs,',',1,1), $rwidth+2);
print_cols("median",  $lwidth, num2str($median,',',1,1),        $rwidth+2);
print_cols("mode",    $lwidth, num2str($mode),                  $rwidth);
print_cols("N50",     $lwidth, num2str($n50),                   $rwidth);

sub print_cols
{
  my ($ltext,$lwidth,$rtext,$rwidth) = @_;
  my $llen = length($ltext);
  my $rlen = length($rtext);
  my ($lpad,$rpad) = ("","");
  if($llen < $lwidth) { $lpad = " "x($lwidth - $llen); }
  if($rlen < $rwidth) {
    $rpad = "."x($rwidth - $rlen);
    substr($rpad,0,1) = substr($rpad,-1,1) = ' ';
  }

  if(!$print_csv) {
    print "[contig_stats.pl] ".$lpad.$ltext.":".$rpad.$rtext."\n";
  } else {
    $rtext =~ s/,//g;
    print "$ltext,$rtext\n";
  }
}

# @_ must be sorted
# If there is a tie for the mode, return the first value
# e.g. 1,2,2,3,4,4 => 2
sub find_mode
{
  if(@_ == 0) { return undef; }
  my ($maxidx,$maxrun,$run) = (0,0,0);
  for(my $i = 1; $i < @_; $i++) {
    if($_[$i] == $_[$i-1]) {
      $run++;
      if($run > $maxrun) { $maxidx = $i; $maxrun = $run; }
    }
    else { $run = 1; }
  }
  return $_[$maxidx];
}

# @_ must be sorted
sub find_median
{
  if(@_ == 0) { return undef; }
  if(@_ % 2 == 0) {
    return ($_[@_ / 2 - 1] + $_[@_ / 2]) / 2;
  } else {
    return $_[int(@_ / 2)];
  }
}

# usage: find_N50($total,@array)
# @array must be sorted
sub find_N50
{
  my $total = shift;
  my ($i, $sum, $half) = (0, 0, $total/2);

  if(@_ == 0 || $half == 0) { return undef; }

  for($i = @_; $i > 0 && $sum < $half; $i--) {
    $sum += $_[$i-1];
  }

  return $_[$i];
}
