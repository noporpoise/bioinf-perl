#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

use FASTNFile;
use GeneticsModule;
use UsefulModule;

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
"Usage: ./sim_reads.pl <OPTIONS> <out_base> <genome.fa>
  options:
    --readlen <len>   [default: $READLEN]
    --mpsize <insert> [default: no matepairs]
    --covg <depth>    [default: $READDEPTH]

  Produces <out_base>.fa or <out_base>.0.fa and <out_base>.1.fa\n";

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
  else { last; }
}

if(@ARGV != 2) { print_usage(); }
my $out_base = shift;
my $genome_file = shift;

my @names = qw(--readlen --mpsize --covg);
my @checks = ($READLEN, $READMPSIZE, $READDEPTH);

for(my $i = 0; $i < @names; $i++) {
  if($checks[$i] !~ /^\d+$/) {
    print_usage("Invalid $names[$i] arg '$checks[$i]'");
  }
}

#
# Load genome
#
my $fastn = open_fastn_file($genome_file);

my $ref = "";
my ($title, $seq);

while((($title,$seq) = $fastn->read_next()) && defined($title))
{
  $seq =~ s/[^ACGT]//gi;
  $ref .= uc($seq);
}

close_fastn_file($fastn);

#
# Sample reads
#
my $genome_len = length($ref);
my $num_reads = $genome_len * $READDEPTH / $READLEN;

my $max_start = $genome_len - ($READMP ? 2*$READLEN+$READMPSIZE : $READLEN) - 1;

if($max_start <= 0) {
  print "Genome too small!\n";
  exit(-1);
}

print "max start: $max_start\n";

if($READMP) {
  my $f0 = "$out_base.0.fa";
  my $f1 = "$out_base.1.fa";
  open(READS0, ">$f0") or die("Cannot write to $f0");
  open(READS1, ">$f1") or die("Cannot write to $f1");
  print "Writing to paired reads to $f0 and $f1\n";
}
else {
  my $f = "$out_base.fa";
  open(READS0, ">$f") or die("Cannot write to $f");
  print "Writing to $f\n";
}

for(my $i = 0; $i < $num_reads; $i++) {
  my $pos = int(rand($max_start));

  my $read = substr($ref, $pos, $READLEN);
  print READS0 ">sample_$i\n$read\n";

  if($READMP) {
    $read = substr($ref, $pos+$READLEN+$READMPSIZE, $READLEN);
    print READS1 ">sample_$i\n" . rev_comp($read) . "\n";
  }
}

close(READS0);
if($READMP) { close(READS0); }
