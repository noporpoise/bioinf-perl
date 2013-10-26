#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule;
use GeneticsModule;
use IndexedString;
use FASTNFile;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }

  print STDERR "Usage: ./sim_substrings.pl <kmer> <approx> <check.fa> <ref1.fa> [ref2.fa] ...
  Count how many of <check.fa> contigs exist in references.
    <kmer> is the seed length used for finding hits
    <approx> is a number >= 0.0 : how much of a length diff is 'approx'\n";
  exit(-1);
}

if(@ARGV < 4) { print_usage(); }

my $kmer_size = shift;
my $approx = shift;
my $check_path = shift;
my @ref_paths = @ARGV;

if($kmer_size !~ /^\d+$/) { print_usage("Invalid kmer size"); }
if($approx !~ /^\d*.?\d*$/) { print_usage("Invalid approx value [must be >= 0.0]"); }

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
my ($num_reads, $num_pass, $num_approx) = (0, 0, 0);
my $fastn = open_fastn_file($check_path);
my ($title,$seq);
while((($title,$seq) = $fastn->read_next()) && defined($title)) {
  $seq =~ s/[^ACGT]//gi;
  $seq = uc($seq);
  my ($idx,$pos) = $search_genomes->find_index($seq);
  if($idx == -1) { ($idx,$pos) = $search_genomes->find_index(rev_comp($seq)); }
  $num_reads++;
  $num_pass += ($idx != -1);

  # Get approx matches
  my $len = length($seq);
  if($idx == -1 && $len > $kmer_size)
  {
    my ($k0, $k1) = (substr($seq, 0, $kmer_size), substr($seq, -$kmer_size));
    if(find_min_dist_diff($k0,$k1,$len) < $len * $approx ||
       find_min_dist_diff(rev_comp($k1),rev_comp($k0),$len) < $len * $approx) {
      $num_approx++;
    }
  }
}
close_fastn_file($fastn);

print "[sim_substrings.pl] Perfect matches: " .
      pretty_fraction($num_pass, $num_reads) . "\n";
print "[sim_substrings.pl] Perfect or approx [".(100*$approx)."%]: " .
      pretty_fraction($num_approx+$num_pass, $num_reads) . "\n";

sub find_min_dist_diff
{
  my ($k0,$k1,$len) = @_;
  my ($kmers0, $kmers1) = map{$search_genomes->{_kmers}->{$_}} ($k0,$k1);
  my @starts = (); my @ends = ();
  add_to_arr($kmers0, \@starts);
  add_to_arr($kmers1, \@ends);
  my $min = LONG_MAX;

  for(my $i = 0; $i < @starts; $i++)
  {
    if(defined($starts[$i]) && defined($ends[$i]) &&
       @{$starts[$i]} > 0 && @{$ends[$i]} > 0)
    {
      $min = find_chrom_min_dist_diff($starts[$i], $ends[$i], $min, $len);
    }
  }

  return $min;
}

sub find_chrom_min_dist_diff
{
  my ($sstarts, $eends, $min, $len) = @_;
  my @starts = sort {$a <=> $b} @$sstarts;
  my @ends = sort {$a <=> $b} @$eends;

  my $e = 0;
  for(my $i = 0; $i < @starts; $i++) {
    while($e < @ends && $ends[$e] < $starts[$i]) { $e++; }
    if($e == @ends) { last; }
    for(my $j = $e; $j < @ends; $j++) {
      my $dist = $ends[$j] - $starts[$i];
      my $diff = abs($len - $dist);
      $min = min($min, $diff);
      if($dist > $len) { last; }
    }
  }

  return $min;
}

sub add_to_arr
{
  my ($list, $arr) = @_;
  for my $k (@$list) {
    my ($idx,$pos) = @$k;
    if(!defined($arr->[$idx])) { $arr->[$idx] = []; }
    push(@{$arr->[$idx]}, $pos);
  }
}
