#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

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

# ------------------------------------------
# Symbol       Meaning     Nucleic Acid
# ------------------------------------------
# A            A           Adenine
# C            C           Cytosine
# G            G           Guanine
# T            T           Thymine
# U            U           Uracil
# M          A or C
# R          A or G
# W          A or T
# S          C or G
# Y          C or T
# K          G or T
# V        A or C or G
# H        A or C or T
# D        A or G or T
# B        C or G or T
# X      G or A or T or C
# N      G or A or T or C

# ------------------------------------------
# Amino acids
# ------------------------------------------
# Ala/A GCU, GCC, GCA, GCG
# Arg/R CGU, CGC, CGA, CGG, AGA, AGG
# Asn/N AAU, AAC
# Asp/D GAU, GAC
# Cys/C UGU, UGC
# Gln/Q CAA, CAG
# Glu/E GAA, GAG
# Gly/G GGU, GGC, GGA, GGG
# His/H CAU, CAC
# Ile/I AUU, AUC, AUA
# Leu/L UUA, UUG, CUU, CUC, CUA, CUG
# Lys/K AAA, AAG
# Met/M AUG
# Phe/F UUU, UUC
# Pro/P CCU, CCC, CCA, CCG
# Ser/S UCU, UCC, UCA, UCG, AGU, AGC
# Thr/T ACU, ACC, ACA, ACG
# Trp/W UGG
# Tyr/Y UAU, UAC
# Val/V GUU, GUC, GUA, GUG
# START AUG
# STOP  UAA, UGA, UAG
