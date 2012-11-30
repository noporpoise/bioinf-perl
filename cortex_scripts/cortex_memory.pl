#!/usr/bin/perl

use strict;
use warnings;

use POSIX qw(ceil);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule; # num2str

#
# Calculate amount of memory required by cortex
#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 29 Jan 2011

sub usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }
  
  print "Usage: ./cortex_memory.pl [options] <mem_height> <mem_width>\n";
  print "  Options: --kmer_size -k <k>   K-mer size of De Bruijn graph [default: 31]\n";
  print "           --colours -c <num>   Number of colours in graph    [default: 1]\n";
  print "           --paths -p <num>     Number of paths in graph      [default: 0]\n";
  exit;
}

if(@ARGV < 2 || @ARGV > 6 || @ARGV % 2 != 0) {
  usage();
}

my $memHeight = $ARGV[@ARGV-2];
my $memWidth = $ARGV[@ARGV-1];

my $kmerSize = 31;
my $numOfColours = 1;
my $numOfPaths = 0;

# Loop through options
for(my $i = 0; $i < @ARGV-2; $i+=2)
{
  if($ARGV[$i] =~ /^-?-k(mer(_?size)?)?$/i)
  {
    $kmerSize = $ARGV[$i+1];

    if($kmerSize !~ /^[0-9]+$/ || $kmerSize < 1 || $kmerSize > 63) {
      usage("invalid kmer");
    }
  }
  elsif($ARGV[$i] =~ /^-?-c(olours?)?$/i)
  {
    $numOfColours = $ARGV[$i+1];
    
    if($numOfColours !~ /^[0-9]+$/ || $numOfColours == 0) {
      usage("invalid number of colours");
    }
  }
  elsif($ARGV[$i] =~ /^-?-p(aths?)?$/i)
  {
    $numOfPaths = $ARGV[$i+1];
    
    if($numOfPaths !~ /^[0-9]+$/) {
      usage("invalid number of colours");
    }

    if($numOfPaths % 8 != 0)
    {
      print STDERR "Warning: paths is usually a multiple of 8 or 0\n";
      # Round up to nearest eight
      $numOfPaths = ceil($numOfPaths / 8) * 8;
    }
  }
  else
  {
    usage();
  }
}

my $numOfPathBytes = $numOfPaths / 8;

if($memHeight !~ /^[0-9]+$/ || $memHeight < 1) {
  usage("invalid memHeight");
}

if($memWidth !~ /^[0-9]+$/ || $memWidth < 1) {
  usage("invalid memHeight");
}

print "width: $memWidth; height: $memHeight; colours: " .
      "$numOfColours; kmer_size: $kmerSize\n";

my $num_of_buckets = 2**$memHeight;
my $num_of_hash_entries = $num_of_buckets * $memWidth;

print "Buckets: " . num2str($num_of_buckets) . "\n";
print "Hashtable entries: " . num2str($num_of_hash_entries) . "\n";

# Round entry size to nearest 8 bytes and multiply by the number of entries
my $element_size = 8*ceil($kmerSize/32) + 5*$numOfColours + 1;

# Add path bytes
$element_size += $numOfColours*$numOfPathBytes;

# Don't need to round up
#$element_size = ceil($element_size / 8) * 8;

my $bytes = $num_of_hash_entries * $element_size;

# ~1GB used for keeping track of bucket usage
$bytes += $num_of_buckets * 2;

print "Memory: " . num2str($bytes) . " bytes (" . mem2str($bytes, 0, 1) . ")\n";

print "Command: cortex_var_" . ($kmerSize > 31 ? "63" : "31") .
      "_c" . $numOfColours . " --kmer_size $kmerSize " .
      "--mem_height $memHeight --mem_width $memWidth\n";
