#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./codons_to_dna.pl <seq>\n";

  exit;
}



my %bases = ('A' => 'GCN', 'R' => 'MGN', 'N' => 'AAY', 'D' => 'GAY',
             'C' => 'TGY', 'Q' => 'CAR', 'E' => 'GAR', 'G' => 'GGN',
             'H' => 'CAY', 'I' => 'ATH', 'L' => 'YUN', 'K' => 'AAR',
             'M' => 'ATG', 'F' => 'TTY', 'P' => 'CCN', 'S' => 'WSN',
             'T' => 'ACN', 'W' => 'TGG', 'Y' => 'TAY', 'V' => 'GTN');

for my $seq (@ARGV)
{
  for(my $i = 0; $i < length($seq); $i++)
  {
    my $codon = substr($seq, $i, 1);

    if(defined(my $dna = $bases{$codon}))
    {
      print $dna;
    }
    else
    {
      print_usage("Invalid codon '$codon'");
    }
  }
  print "\n";
}
