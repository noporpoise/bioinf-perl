#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

use FASTNFile;

# config
#use constant {NUM_SNPS => 10000, NUM_INDELS => 2000,
use constant {NUM_SNPS => 100000, NUM_INDELS => 10000, NUM_INV => 1000,
              READDEPTH => 10, READLEN => 100, READMP => 1, READMPSIZE => 250};

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./sim_mutations.pl <kmer_size> <in.fa> <outdir>\n";
  print STDERR "  Creates genome.0.fa, genome.1.fa,\n";
  print STDERR "          reads.0.fa, reads.1.fa,\n";
  print STDERR "          truth.vcf\n";

  exit(-1);
}

if(@ARGV != 3) { print_usage(); }
my $kmer_size = shift;
my $genome_file = shift;
my $outdir = shift;

if(!(-e $outdir)) {
  make_path($outdir) or die("Cannot make dir $outdir");
} elsif(!(-d $outdir)) {
  print_usage("Cannot create output directory - file with same name exists");
}

if($kmer_size !~ /^\d+$/) { print_usage(); }

my $fastn = open_fastn_file($genome_file);

my $ref = "";

my ($title,$seq);

while((($title,$seq) = $fastn->read_next()) && defined($title))
{
  $seq =~ s/[^ACGT]//gi;
  $ref .= uc($seq);
}

close_fastn_file($fastn);

my @genomes = ($ref, $ref);
my $reflen = length($ref);

print "Genome size: $reflen\n";

if($reflen <= NUM_SNPS+NUM_INDELS) {
  print "Warning: genome is smaller than number of snps+indels " .
        "(".NUM_SNPS."+".NUM_INDELS.")\n";
}

# In the mask: . = ref, N = variant
my $mask = '.' x $reflen;

my $num_of_snps = 0;
my $num_of_indels = 0;
my $num_of_invs = 0;

# Variants are non-overlapping
# Generate SNPs
for(my $i = 0; $i < NUM_SNPS; $i++)
{
  my $pos = int(rand($reflen));

  if(substr($mask, $pos, 1) eq ".")
  {
    # Pick a random nonref base (ACGT)
    my $ref_base = substr($genomes[0], $pos, 1);
    my @bases = grep {$_ ne $ref_base} qw(A C G T);
    my $snp = $bases[int(rand(3))];

    substr($genomes[1], $pos, 1) = $snp;
    substr($mask, $pos, 1) = 'N';
    $num_of_snps++;
  }
}

# Generate indels with geometric size distribution
for(my $i = 0; $i < NUM_INDELS; $i++)
{
  my $pos = int(rand($reflen));

  if(substr($mask, $pos, 1) eq '.')
  {
    my $end = $pos;
    while($end+1 < $reflen && rand() < 0.5 && substr($mask, $end+1, 1) eq '.'){$end++;}
    my $len = $end-$pos+1;

    substr($genomes[int(rand(2))], $pos, $len) = '-'x$len;
    substr($mask, $pos, $len) = 'N'x$len;
    $num_of_indels++;
  }
}

# Generate inversions with geometric size distribution
for(my $i = 0; $i < NUM_INV; $i++)
{
  my $pos = int(rand($reflen));

  if(substr($mask, $pos, 1) eq '.')
  {
    my $end = $pos;
    while($end+1 < $reflen && rand() < 0.5 && substr($mask, $end+1, 1) eq '.'){$end++;}
    my $len = $end-$pos+1;

    my $g0 = int(rand(2)); # 0 or 1
    my $g1 = 1 - $g0; # 1 or 0
    substr($genomes[$g0], $pos, $len) = rev_cmp(substr($genomes[$g1], $pos, $len));
    substr($mask, $pos, $len) = 'N'x$len;
    $num_of_invs++;
  }
}


# Print truth VCF
open(VCF, ">$outdir/truth.vcf") or die("Cannot open $outdir/truth.vcf");

my $left = ('.'x$kmer_size).'N';
my $right = 'N'.('.'x$kmer_size);
my $m1 = 0;
my $m2 = 0;

for(my $i = 0;
    ($m1 = index($mask, $left, $m2)) != -1 && ($m2 = index($mask, $right, $m1)) != -1;
    $i++)
{
  my $lflank = substr($ref, $m1, $kmer_size);
  my $rflank = substr($ref, $m2+1, $kmer_size);
  my $start = $m1+$kmer_size;
  my $len = $m2-$start+1;
  my @alleles = map {substr($genomes[$_], $start, $len)} 0..1;
  map {$alleles[$_] =~ s/\-//g} 0..1;

  if(defined(first {length($_) != 1} @alleles))
  {
    # Not a SNP
    my $lastbase = substr($lflank,-1);
    map {$alleles[$_] = $lastbase.$alleles[$_]} 0..1;
  }

  print VCF ".\t0\tvar$i\tN\t" . join(',', @alleles) . "\t" .
            "LF=$lflank;RF=$rflank;BN=2\t.\n";
}

close(VCF);

# Print genomes to genome.0.fa and genome.1.fa
open(GENOM, ">$outdir/genome.0.fa") or die("Cannot write to $outdir/genome.0.fa");
print GENOM ">genomes[0]\n$genomes[0]\n";
close(GENOM);
open(GENOM, ">$outdir/genome.1.fa") or die("Cannot write to $outdir/genome.1.fa");
print GENOM ">genomes[1]\n$genomes[1]\n";
close(GENOM);

# Dump ref
open(GENOM, ">$outdir/ref.fa") or die("Cannot write to $outdir/ref.fa");
print GENOM ">ref\n$ref\n";
close(GENOM);

# Dump mask
open(GENOM, ">$outdir/mask.fa") or die("Cannot write to $outdir/mask.fa");
print GENOM ">mask\n$mask\n";
close(GENOM);

$genomes[0] =~ s/-//g;
$genomes[1] =~ s/-//g;

# Sample reads
my $mean_genome_len = (length($genomes[0]) + length($genomes[1])) / 2;
my $num_reads = $mean_genome_len * READDEPTH / READLEN;

for(my $i = 0; $i < 2; $i++) {
  if(READMP) {
    open(READS0, ">$outdir/reads$i.0.fa") or die("Cannot write to $outdir/reads$i.0.fa");
    open(READS1, ">$outdir/reads$i.1.fa") or die("Cannot write to $outdir/reads$i.1.fa");
  }
  else {
    open(READS0, ">$outdir/reads.$i.fa") or die("Cannot write to $outdir/reads.$i.fa");
  }

  my $genlen = length($genomes[$i]) - (READMP ? 2*READLEN+READMPSIZE : READLEN);

  for(my $j = 0; $j < $num_reads; $j++) {
    my $pos = int(rand($genlen));
    print READS0 ">sample".$i."_$j\n".substr($genomes[$i], $pos, READLEN)."\n";

    if(READMP) {
      print READS1 ">sample".$i."_$j\n" .
                   substr($genomes[$i], $pos+READLEN+READMPSIZE, READLEN) . "\n";
    }
  }

  close(READS0);
  if(READMP) { close(READS0); }
}

print " snps: $num_of_snps / ".NUM_SNPS." generated\n";
print " indels: $num_of_indels / ".NUM_INDELS." generated\n";
print " inversions: $num_of_invs / ".NUM_INV." generated\n";
