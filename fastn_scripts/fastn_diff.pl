#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use GeneticsModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./fastn_diff.pl [--no-rev-cmp] <file1> <file2>
Print sequences in <file2> not in <file1>.

  --no-rev-cmp      Do not compare including reverse complement
  --case-sensitive  Do case sensitive comparision
\n";

  exit;
}

my $case_sensitive = 0;
my $no_revcmp = 0;

while(@ARGV > 2) {
  my $arg = shift(@ARGV);
  if($arg eq "--no-rev-cmp") { $no_revcmp = 1; }
  elsif($arg eq "--case-sensitive") { $case_sensitive = 1; }
  else { print_usage("Bad option: $arg"); }
}

if(@ARGV != 2) { print_usage(); }

my ($path1,$path2) = @ARGV;

my $fastn1 = open_fastn_file($path1);
my $fastn2 = open_fastn_file($path2);

my %seqhash = ();
my ($title,$seq,$qual);

# dna_rev_comp_group() is kmer-key

# Load first file sequences
while((($title, $seq) = $fastn1->read_next()) && defined($title))
{
  if(!$case_sensitive) { $seq = uc($seq); }
  if(!$no_revcmp) { $seq = dna_rev_comp_group($seq); }
  $seqhash{$seq} = 1;
}

# Load second file
while((($title, $seq, $qual) = $fastn2->read_next()) && defined($title))
{
  my $cpy = $seq;
  if(!$case_sensitive) { $cpy = uc($cpy); }
  if(!$no_revcmp) { $cpy = dna_rev_comp_group($cpy); }
  if(!defined($seqhash{$cpy})) {
    $fastn2->print_entry($title,$seq,$qual);
  }
}
