#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first sum shuffle);
use File::Path qw(make_path);

use FASTNFile;
use GeneticsModule;
use UsefulModule;

# config defaults
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
  
  print STDERR "" .
"Usage: ./sim_mutations.pl [options] <outdir> <num_samples> [in.fa|in.fq.gz]\n" .
"  options: --snps <num> --indels <num> --invs <num> --invlen <len>\n" .
"  Creates: ref.fa, genome1.fa .. genomeN.fa, mask1.fa .. maskN.fa
    where N is <num_samples>\n";

  exit(-1);
}

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^--snps$/i) {
    shift;
    $NUM_SNPS = shift;
  }
  if($ARGV[0] =~ /^--indels$/i) {
    shift;
    $NUM_INDELS = shift;
  }
  if($ARGV[0] =~ /^--invs$/i) {
    shift;
    $NUM_INV = shift;
  }
  if($ARGV[0] =~ /^--invlen$/i) {
    shift;
    $INVLEN = shift;
  }
  else { last; }
}

if(@ARGV < 2 || @ARGV > 3) { print_usage(); }
my $outdir = shift;
my $num_of_samples = shift;
my $genome_file = shift;

if(!defined($genome_file)) { $genome_file = "-"; }

if($num_of_samples !~ /^\d+$/ || $num_of_samples == 0) {
  print_usage("<num_samples> must be a +ve int");
}

my @names = qw(--snps --indels --inv --invlen);
my @checks = ($NUM_SNPS, $NUM_INDELS, $NUM_INV, $INVLEN);

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

#
# Load genome
#
my $fastn = open_fastn_file($genome_file);
my $ref = "";
my ($title,$seq);

while((($title,$seq) = $fastn->read_next()) && defined($title))
{
  $seq =~ s/[^ACGT]//gi;
  $ref .= uc($seq);
}

close_fastn_file($fastn);


my $reflen = length($ref);
print "Genome size: ".num2str($reflen)."\n";

if($READMP && $reflen < 2*$READLEN+$READMPSIZE) {
  print "Error: genome is smaller than `read + mpinsertsize + read`\n";
  exit(-1);
}
if($reflen <= $NUM_SNPS+$NUM_INDELS+$NUM_INV) {
  print "Warning: genome is smaller than number of snps+indels (" .
        join(" + ", num2str($NUM_SNPS), num2str($NUM_INDELS), num2str($NUM_INV)).
        ")\n";
}

# In the mask: . = ref, N = variant
my @genomes = ($ref) x $num_of_samples;
my @masks = ('.' x $reflen) x $num_of_samples;
my $mask = '.' x $reflen;

my $num_of_snps = 0;
my $num_of_ins = 0;
my $num_of_del = 0;
my $num_of_invs = 0;

# Generate likelihood of N samples having mutations
my @prob;
if($num_of_samples == 1) { $prob[0] = 1; }
else {
  @prob = (0, map {1 / $_} 1..($num_of_samples-1));
  # Normalise
  my $prob_sum = sum(@prob);
  @prob = map {$_ / $prob_sum} @prob;
  # Convert to cdf
  for(my $i = 1; $i < $num_of_samples; $i++) { $prob[$i] += $prob[$i-1]; }
}

# print "".join(', ', map {$_.":".$prob[$_]} 0..($num_of_samples-1))."\n";

my @sampleids = 0..($num_of_samples-1);

# Variants are non-overlapping
# Generate SNPs
for(my $i = 0; $i < $NUM_SNPS; $i++)
{
  my $pos = int(rand($reflen-1)+1);

  if(substr($mask, $pos, 1) eq ".")
  {
    # Pick a random nonref base (ACGT)
    my $ref_base = substr($ref, $pos, 1);
    my @bases = grep {$_ ne $ref_base} qw(A C G T);
    my $snp = $bases[int(rand(3))];

    my @snp_samples = (0);
    if($num_of_samples > 1) {
      # Number of affected individuals (integer {1..N-1} distib 1/N)
      my $p = rand();
      my $N = first {$p < $prob[$_]} 1..($num_of_samples-1);
      my @samples = shuffle(@sampleids);
      @snp_samples = @samples[0..($N-1)];
    }

    for my $s (@snp_samples)
    {
      substr($genomes[$s], $pos, 1) = $snp;
      substr($masks[$s], $pos, 1) = 'S';
    }

    # Pick one to make lowercase
    substr($masks[$snp_samples[0]], $pos, 1) = 's';

    substr($mask, $pos, 1) = 'S';
    $num_of_snps++;
  }
}

# Generate indels with geometric size distribution
for(my $i = 0; $i < $NUM_INDELS; $i++)
{
  my $pos = int(rand($reflen-1)+1);

  if(substr($mask, $pos, 1) eq '.')
  {
    # generate length (geometric distribution)
    my $end = $pos;
    while($end+1 < $reflen && rand() < 0.5 && substr($mask, $end+1, 1) eq '.'){$end++;}
    my $len = $end-$pos+1;


    my $is_ins = 0;
    my @del_samples = (0);
    my @ins_samples = ();

    if($num_of_samples > 1) {
      # insertion or deletion
      $is_ins = rand() < 0.5;

      # Number of affected individuals (integer {1..N-1} distib 1/N)
      my $p = rand();
      my $N = first {$p < $prob[$_]} 1..($num_of_samples-1);

      # make N the number of people we're 'deleting' from
      if($is_ins) { $N = $num_of_samples - $N; }

      # Pick random individuals
      my @samples = shuffle(@sampleids);
      @del_samples = @samples[0..($N-1)];
      @ins_samples = @samples[$N..($num_of_samples-1)];
    }

    # delete/insert from selected individuals
    for my $s (@del_samples) {
      substr($genomes[$s], $pos, $len) = '-'x$len;
    }

    if($is_ins) {
      for my $s (@ins_samples) {
        substr($masks[$s], $pos, $len) = 'I'x$len;
      }
      # Pick one to be lower case at start
      substr($masks[$ins_samples[0]], $pos, 1) = 'i';
      $num_of_ins++;
    }
    else {
      for my $s (@del_samples) {
        substr($masks[$s], $pos, $len) = 'D'x$len;
      }
      # Pick one to be lower case at start
      substr($masks[$del_samples[0]], $pos, 1) = 'd';
      $num_of_del++;
    }

    substr($mask, $pos, $len) = 'N'x$len;
  }
}

# Generate inversions with uniform size distribution Uniform(5,INVLEN) bp long
for(my $i = 0; $i < $NUM_INV; $i++)
{
  my $pos = int(rand($reflen-1)+1);
  my $len = int(rand($INVLEN-5))+5;

  if($pos + $len > $reflen) { $len = $reflen - $pos; }

  if(substr($mask, $pos, $len) =~ /^\.+$/)
  {
    my @inv_samples = (0);

    if($num_of_samples > 0)
    {
      # Number of affected individuals (integer {1..N-1} distib 1/N)
      my $p = rand();
      my $N = first {$p < $prob[$_]} 1..($num_of_samples-1);

      my @samples = shuffle(@sampleids);
      @inv_samples = @samples[0..($N-1)];
    }

    my $inv = rev_comp(substr($ref, $pos, $len));
    for my $s (@inv_samples)
    {
      substr($genomes[$s], $pos, $len) = $inv;
      substr($masks[$s], $pos, $len) = 'V'x$len;
    }

    # Pick one to be lower case at start
    substr($masks[$inv_samples[0]], $pos, 1) = 'v';

    substr($mask, $pos, $len) = 'V'x$len;
    $num_of_invs++;
  }
}


# Print genomes to genome.0.fa and genome.1.fa
my $f;
for(my $i = 0; $i < $num_of_samples; $i++)
{
  $f = "$outdir/genome$i.fa";
  open(GENOM, ">$f") or die("Cannot write to $f");
  print GENOM ">genome$i\n".$genomes[$i]."\n";
  close(GENOM);

  $f = "$outdir/mask$i.fa";
  open(MASK, ">$f") or die("Cannot write to $f");
  print MASK ">mask$i\n".$masks[$i]."\n";
  close(MASK);
}

print " snps: ".pretty_fraction($num_of_snps, $NUM_SNPS)." generated\n";
print " insertions: ".pretty_fraction($num_of_ins, $NUM_INDELS/2)." generated\n";
print " deletions: ".pretty_fraction($num_of_del, $NUM_INDELS/2)." generated\n";
print " inversions: ".pretty_fraction($num_of_invs, $NUM_INV)." generated\n";
