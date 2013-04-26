#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

my $csvsep = "\t";

# Chimps
#my $read_depth = 6;
#my $num_samples = 10;
#my $read_length_bp = 51;
#my $kmer_size = 31;

# Simon's Mice
#my $read_depth = 10;
#my $num_samples = 20;
#my $read_length_bp = 100;
#my $kmer_size = 31;

#my $ploidy = 2;
#my $genome_complexity = 0.75;
#my $epsilon = 0.005;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./calculate_power.pl <max_indel> <depth> <ploidy> " .
               "<samples> <readlen> <kmer> <complexity> <epsilon>
  Prints discovery power matrix (indel_size x allele freq)\n";
  exit;
}

if(@ARGV != 8)
{
  print_usage();
}

my $max_indel = shift;
my $read_depth = shift;
my $ploidy = shift;
my $num_samples = shift;
my $read_length_bp = shift;
my $kmer_size = shift;
my $genome_complexity = shift;
my $epsilon = shift;

my @int_values = ($max_indel, $ploidy, $num_samples, $read_length_bp, $kmer_size);
my @int_names = qw(max_indel ploidy num_samples read_length_bp kmer_size);

for(my $i = 0; $i < @int_names; $i++)
{
  if($int_values[$i] !~ /^\d+$/)
  {
    print_usage("$int_names[$i] must be a positive integer ('$int_values[$i]' is invalid)");
  }
}

my @double_values = ($read_depth, $genome_complexity, $epsilon);
my @double_names = qw(read_depth genome_complexity epsilon);

for(my $i = 0; $i < @double_names; $i++)
{
  if($double_values[$i] !~ /^\d+(\.\d+)?$/)
  {
    print_usage("$double_names[$i] must be a number ('$double_values[$i]' is invalid)");
  }
}

if($genome_complexity > 1)
{
  print_usage("genome_complexity cannot be > 1");
}

if($epsilon > 1)
{
  print_usage("epsilon cannot be > 1");
}

# Covg depth is coverage per chrom
my $covg_depth = $read_depth / $ploidy;
my $read_kmer_length = $read_length_bp - $kmer_size + 1; # = L

#my @freqs = (0.05, (map {$_/10} 1..18), 1.95);

my @freqs = reverse(0.025, 0.05, (map {$_/10} 1..9), 0.95, 0.975);
my @indels = (1,5,10,15,20,25,30);

my $col_headings = join($csvsep, (map {"-$_"} reverse(@indels)), @indels);

#my @freqs = (0.05);

open(OUT, ">power.csv") or die("Can't open power.csv");

#non_ref_freq
print "0".$csvsep.$col_headings."\n";
print OUT "0".$csvsep.$col_headings."\n";

for my $freq (@freqs)
{
  my $line = $freq;

  # Deletions from the reference
  for my $i (reverse(@indels))
  {
    # $freq 31 kmers vs [1-$freq] 31+$i
    my $power = get_power(31, $freq, 31+$i, 1-$freq);

    $line .= $csvsep . sprintf("%.3f", $power);
  }

  # Insertions to the reference
  for my $i (@indels)
  {
    my $power = get_power(31+$i, $freq, 31, 1-$freq);

    $line .= $csvsep . sprintf("%.3f", $power);
  }

  $line .= "\n";

  print $line;
  print OUT $line;
}

close(OUT);

print "0.05 means one of 10 individuals is Aa, rest are AA\n";

sub get_power
{
  my ($len_a, $freq_a, $len_b, $freq_b) = @_;

  return $genome_complexity * (1- $kmer_size * $epsilon) *
         allele_covg_prob($len_a, $freq_a) * allele_covg_prob($len_b, $freq_b);
}

# $allele length is var size + k bp either side
sub allele_covg_prob
{
  my ($allele_length_bp, $freq) = @_;

  # Round allele length the nearest 10bp to get the number of supernodes that
  # are as long or longer than required
  #my $supernode_count = $total_nodes_gt[round_int($allele_length_bp+30, 10) / 10];

  # Calculate power
  my $lambda = (2 * $num_samples * $freq * $covg_depth) / $read_length_bp;

  # t = allele length
  #return ($supernode_count / $total_count) * 
  my $power = ((1 - exp(-$lambda * $read_kmer_length)) ** 2) *
              exp(-$lambda * $allele_length_bp * exp(-$lambda * $read_kmer_length));

  # print "((1 - exp(-$lambda * $read_kmer_length)) ^ 2) * " .
  #       "exp(-$lambda * $allele_length_bp * " .
  #       "exp(-$lambda * $read_kmer_length)) = $power\n";

  #print "$allele_length_bp bp; $freq freq; lambda: $lambda; power: $power\n";

  return $power;
}
