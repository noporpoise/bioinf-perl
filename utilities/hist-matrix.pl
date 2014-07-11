#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(reduce);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./hist-matrix.pl <dim> <file>
  Converts 3-column input into 3-dimensional matrix of counts
\n";
  exit(-1);
}

if(@ARGV != 2) { print_usage(); }

my $dim = shift(@ARGV);
my $path = shift(@ARGV);

if($dim !~ /^\d+$/) { print_usage(); }

my $dim2 = $dim *$dim;
my $dim3 = $dim2*$dim;
my $idx;
my @arr = (0)x$dim3;

my $fh;
my $line;
if(-p STDIN) { open($fh, "<&=STDIN") or die("Cannot read pipe"); }
else         { open($fh, $path) or die("Cannot open file: $path"); }

while(defined($line = <$fh>))
{
  my @cols = map {$_ >= $dim ? $dim-1 : $_-1} split(/[\s,]/, $line);
  # Matrices have backwards indices so we gotta do backwards indices
  my $idx = $cols[2]*$dim2 + $cols[0]*$dim + $cols[1];
  $arr[$idx]++;
}

close($fh);

$idx = 0;
for(my $i = 0; $i < $dim; $i++) {
  for(my $j = 0; $j < $dim; $j++) {
    print "".join(',', @arr[$idx..($idx+$dim-1)])."\n";
    $idx += $dim;
  }
}
