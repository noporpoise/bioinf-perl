#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

use FASTNFile;
use GeneticsModule;
use UsefulModule;

# config
# use constant {NUM_SNPS => 100000, NUM_INDELS => 10000, NUM_INV => 1000,
              # READDEPTH => 10, READLEN => 100, READMP => 1, READMPSIZE => 450};

my $NUM_SNPS = 100000;
my $NUM_INDELS = 10000;
my $NUM_INV = 1000;
my $INVLEN = 500;
my $READDEPTH = 10;
my $READLEN = 100;
my $READMP = 0;
my $READMPSIZE = 450;

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./sim_mutations.pl [options] <kmer_size> <in.fa> <outdir>\n";
  print STDERR "  options: --readlen <len> --mpsize <insert> --covg <depth>\n";
  print STDERR "           --snps <num> --indels <num> --inv <num> --invlen <len>\n";
  print STDERR "  Creates genomeA.fa, genomeB.fa, ref.fa, mask.fa, truth.vcf\n";
  print STDERR "          readsA.fa, readsB.fa OR if --mpsize given:\n";
  print STDERR "          readsA.0.fa, readsA.1.fa, readsB.0.fa readsB.1.fa\n";

  exit(-1);
}

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^--readlen$/i) {
    shift;
    $READLEN = shift;
  }
  if($ARGV[0] =~ /^--mpsize$/i) {
    shift;
    $READMP = 1;
    $READMPSIZE = shift;
  }
  if($ARGV[0] =~ /^--covg$/i) {
    shift;
    $READDEPTH = shift;
  }
  if($ARGV[0] =~ /^--snps$/i) {
    shift;
    $NUM_SNPS = shift;
  }
  if($ARGV[0] =~ /^--indels$/i) {
    shift;
    $NUM_INDELS = shift;
  }
  if($ARGV[0] =~ /^--inv$/i) {
    shift;
    $NUM_INV = shift;
  }
  if($ARGV[0] =~ /^--invlen$/i) {
    shift;
    $INVLEN = shift;
  }
  else { last; }
}

if(@ARGV != 3) { print_usage(); }
my $kmer_size = shift;
my $genome_file = shift;
my $outdir = shift;

my @names = qw(--readlen --mpsize --covg --snps --indels --inv --invlen);
my @checks = ($READLEN, $READMPSIZE, $READDEPTH,
              $NUM_SNPS, $NUM_INDELS, $NUM_INV, $INVLEN);

for(my $i = 0; $i < @names; $i++) {
  if($checks[$i] !~ /^\d+$/) {
    print_usage("Invalid $names[$i] arg '$checks[$i]'");
  }
}

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

if($READMP && $reflen < 2*$READLEN+$READMPSIZE) {
  print "Error: genome is smaller than `read + mpinsertsize + read`\n";
  exit(-1);
}
if($reflen <= $NUM_SNPS+$NUM_INDELS+$NUM_INV) {
  print "Warning: genome is smaller than number of snps+indels " .
        "(".num2str($NUM_SNPS)."+".num2str($NUM_INDELS)."+".num2str($NUM_INV).")\n";
}

# In the mask: . = ref, N = variant
my $mask = '.' x $reflen;

my $num_of_snps = 0;
my $num_of_indels = 0;
my $num_of_invs = 0;

# Variants are non-overlapping
# Generate SNPs
for(my $i = 0; $i < $NUM_SNPS; $i++)
{
  my $pos = int(rand($reflen));

  if(substr($mask, $pos, 1) eq ".")
  {
    # Pick a random nonref base (ACGT)
    my $ref_base = substr($genomes[0], $pos, 1);
    my @bases = grep {$_ ne $ref_base} qw(A C G T);
    my $snp = $bases[int(rand(3))];

    substr($genomes[1], $pos, 1) = $snp;
    substr($mask, $pos, 1) = 'S';
    $num_of_snps++;
  }
}

# Generate indels with geometric size distribution
for(my $i = 0; $i < $NUM_INDELS; $i++)
{
  my $pos = int(rand($reflen));

  if(substr($mask, $pos, 1) eq '.')
  {
    my $end = $pos;
    while($end+1 < $reflen && rand() < 0.5 && substr($mask, $end+1, 1) eq '.'){$end++;}
    my $len = $end-$pos+1;

    substr($genomes[int(rand(2))], $pos, $len) = '-'x$len;
    substr($mask, $pos, $len) = 'I'x$len;
    $num_of_indels++;
  }
}

# Generate inversions with uniform size distribution 5 -> INVLEN bp
for(my $i = 0; $i < $NUM_INV; $i++)
{
  my $pos = int(rand($reflen));
  my $len = int(rand($INVLEN-5))+5;

  if(substr($mask, $pos, $len) =~ /^\.+$/)
  {
    my $g0 = int(rand(2)); # 0 or 1
    my $g1 = 1 - $g0; # 1 or 0
    substr($genomes[$g0], $pos, $len) = rev_comp(substr($genomes[$g1], $pos, $len));
    substr($mask, $pos, $len) = 'V'x$len;
    $num_of_invs++;
  }
}


# Print truth VCF
open(VCF, ">$outdir/truth.vcf") or die("Cannot open $outdir/truth.vcf");

my $left = ('.'x$kmer_size).'N';
my $right = 'N'.('.'x$kmer_size);
my $m1 = 0;
my $m2 = 0;

for(my $i = 0; $mask =~ /(\.{$kmer_size})[^\.]+(\.{$kmer_size})/g; $i++)
{
  # $-[1] is pos of left flank, $-[2] is pos of right flank
  my $start = $-[0]+$kmer_size;
  my $len = $-[2] - $-[1] - $kmer_size;

  my $lflank = substr($ref, $-[1], $kmer_size);
  my $rflank = substr($ref, $-[2], $kmer_size);

  my $m = substr($mask, $start, $len);
  my @alleles = map {substr($genomes[$_], $start, $len)} 0..1;

  map {$alleles[$_] =~ s/\-//g} 0..1;

  if(defined(first {length($_) != 1} @alleles))
  {
    # Not a SNP
    map {$alleles[$_] = 'N'.$alleles[$_]} 0..1;
  }

  print VCF ".\t0\tvar$i\tN\t" . join(',', @alleles) . "\t" .
            "LF=$lflank;RF=$rflank";

  if($m =~ /S/) { print VCF ";SNP"; }
  if($m =~ /I/) { print VCF ";INDEL"; }
  if($m =~ /V/) { print VCF ";INV"; }

  print VCF "\t.\n";
}

close(VCF);

# Print genomes to genome.0.fa and genome.1.fa
open(GENOM, ">$outdir/genomeA.fa") or die("Cannot write to $outdir/genomeA.fa");
print GENOM ">genomes[0]\n$genomes[0]\n";
close(GENOM);
open(GENOM, ">$outdir/genomeB.fa") or die("Cannot write to $outdir/genomeB.fa");
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
my $num_reads = $mean_genome_len * $READDEPTH / $READLEN;

for(my $i = 0; $i < 2; $i++) {
  my $indv = ($i == 0 ? 'A' : 'B');
  if($READMP) {
    open(READS0, ">$outdir/reads$indv.0.fa") or die("Cannot write to $outdir/reads$indv.0.fa");
    open(READS1, ">$outdir/reads$indv.1.fa") or die("Cannot write to $outdir/reads$indv.1.fa");
  }
  else {
    open(READS0, ">$outdir/reads$indv.fa") or die("Cannot write to $outdir/reads$indv.fa");
  }

  my $genlen = length($genomes[$i]) - ($READMP ? 2*$READLEN+$READMPSIZE : $READLEN);

  for(my $j = 0; $j < $num_reads; $j++) {
    my $pos = int(rand($genlen));

    my $read = substr($genomes[$i], $pos, $READLEN);
    print READS0 ">sample".$i."_$j\n$read\n";

    if($READMP) {
      $read = substr($genomes[$i], $pos+$READLEN+$READMPSIZE, $READLEN);
      print READS1 ">sample".$i."_$j\n" . rev_comp($read) . "\n";
    }
  }

  close(READS0);
  if($READMP) { close(READS0); }
}

print " snps: ".num2str($num_of_snps)." / ".num2str($NUM_SNPS)." generated\n";
print " indels: ".num2str($num_of_indels)." / ".num2str($NUM_INDELS)." generated\n";
print " inversions: ".num2str($num_of_invs)." / ".num2str($NUM_INV)." generated\n";
