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
  my $idx = $cols[1]*$dim2 + $cols[0]*$dim + $cols[2]; # for R matrix
  # my $idx = $cols[2]*$dim2 + $cols[0]*$dim + $cols[1];
  $arr[$idx]++;
}

close($fh);

# Input is: [num_kmers] [num_juncs] [score]

# Output is 3-dim matrix:
# block/depth is [num_juncs]
# row is [num_kmers]
# column is [score]

# column, row and blocks start at 1

$idx = 0;
for(my $i = 0; $i < $dim; $i++) {
  for(my $j = 0; $j < $dim; $j++) {
    print "".join(',', @arr[$idx..($idx+$dim-1)])."\n";
    $idx += $dim;
  }
}

# Load into R with
#   x=read.csv(file='mat.txt',header=F)
#   m=as.matrix(x)
#   dim(m)=c(3,3,3)
#
# Address with:
#    m[num_kmers,num_juncs,score] = count
#

# Scratch space
# z=abind(x[1:3,1:3],x[4:6,1:3],x[7:9,1:3],along=3)
