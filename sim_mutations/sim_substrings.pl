#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule;
use IndexedString;
use FASTNFile;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "Usage: ./sim_substrings.pl <check.fa> <ref1.fa> [ref2.fa] ...
  Count how many of <check.fa> contigs exist in references.\n";
  exit(-1);
}

if(@ARGV < 2) { print_usage(); }

my $kmer_size = 31;

my $check_path = shift;
my @ref_paths = @ARGV;

my $search_genomes = new IndexedString($kmer_size);

for my $ref_path (@ref_paths)
{
  my $fastn = open_fastn_file($ref_path);
  my ($title,$seq);
  while((($title,$seq) = $fastn->read_next()) && defined($title)) {
    $seq =~ s/[^ACGT]//gi;
    $seq = uc($seq);
    $search_genomes->add_strings($seq);
  }
  close_fastn_file($fastn);
}

# Read input file
my ($num_reads, $num_pass) = (0, 0);
my $fastn = open_fastn_file($check_path);
my ($title,$seq);
while((($title,$seq) = $fastn->read_next()) && defined($title)) {
  $seq =~ s/[^ACGT]//gi;
  $seq = uc($seq);
  my ($idx,$pos) = $search_genomes->find_index($seq);
  $num_reads++;
  $num_pass += ($idx != -1);
}
close_fastn_file($fastn);

print "Perfect matches: ".pretty_fraction($num_pass, $num_reads)."\n";
