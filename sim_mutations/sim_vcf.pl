#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

use FASTNFile;
use GeneticsModule;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./sim_vcf.pl <kmer_size> [<sample.fa> <sample.mask> ..]\n";

  exit(-1);
}

if(@ARGV < 3) { print_usage(); }

my $kmer_size = shift;

if($kmer_size !~ /^\d+$/) { print_usage(); }
if(scalar(@ARGV) % 2 != 0) { print_usage(); }

# All masks and chromosomes should be the same length

my @genomes = ();
my @masks = ();

#
# Load and merge mask files
#
my $num_of_samples = int(@ARGV / 2);

for(my $i = 0; $i < $num_of_samples; $i++)
{
  my $genome_file = shift;
  my $mask_file = shift;

  my $fastn;
  my ($title,$seq);

  $fastn = open_fastn_file($genome_file);
  ($title,$seq) = $fastn->read_next();
  close_fastn_file($fastn);
  if(!defined($title)) { die("Empty file: $genome_file\n"); }
  push(@genomes, uc($seq));

  $fastn = open_fastn_file($mask_file);
  ($title,$seq) = $fastn->read_next();
  close_fastn_file($fastn);
  if(!defined($title)) { die("Empty file: $mask_file\n"); }
  push(@masks, $seq);
}

#
# Merge masks
#
my $mask = "." x length($masks[0]);

for(my $i = 0; $i < $num_of_samples; $i++)
{
  while($masks[$i] =~ /([^\.]+)/)
  {
    my $start = $-[0];
    my $end = $+[0];
    substr($mask, $start, $end-$start+1) = $1;
  }
}

#
# Remove '-' that are shared by all individuals
#
my $pos = 0;
while(($pos = index($mask, $pos)) > -1)
{
  my $len = 0;
  while(defined(first {$_ ne '-'} map {substr($_,$pos,$len+1)} @masks)) {
    $len++;
  }
  if($len > 0) {
    for(my $i = 0; $i < $num_of_samples; $i++) {
      substr($masks[$i],$pos,$len) = 0;
    }
    substr($mask,$pos,$len) = 0;
  }
}

#
# Generate truth VCF
#
for(my $i = 0; $mask =~ /\.{$kmer_size}([^\.]+)(\.{$kmer_size})/g; $i++)
{
  # $-[1] is pos of left flank, $-[2] is pos of right flank
  my $start = $-[0];
  my $len = length($1);

  # Get alleles
  my @alleles = map {$_ =~ s/\-//g}
                map {substr($genomes[$_], $start, $len)}
                0..$num_of_samples;

  # Remove duplicates
  my %ah = ();
  @ah{@alleles} = 1;
  @alleles = keys(%ah);

  if(@alleles > 1)
  {
    # Trim matching bp off ends
    my ($lend,$rend) = trim_alleles(@alleles);
    $start += $lend;
    $len -= $lend + $rend;

    my $lflank = substr($genomes[0], $start-$kmer_size, $kmer_size);
    my $rflank = substr($genomes[0], $start+$len, $kmer_size);

    if(defined(first {length($_) != 1} @alleles))
    {
      # Not a SNP
      map {$alleles[$_] = 'N'.$alleles[$_]} 0..1;
    }

    my $info = "LF=$lflank;RF=$rflank";

    # Get annotations from the mask
    my $m = substr($mask, $start, $len);
    if($m =~ /S/) { $info .= ";SNP"; }
    if($m =~ /I/i) { $info .= ";INDEL"; }
    if($m =~ /V/) { $info .= ";INV"; }

    print ".\t0\tvar$i\tN\t" . join(',', @alleles) . "\t$info\t.\n";
  }
}

# Get number of matching bases on left and right
sub trim_alleles
{
  my @a = @_;
  my ($l, $r) = (0, 0);

  for(; $l < length($_[0]); $l++) {
    my @bp = map {substr($_,$l,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  my $remaining = lenth($_[0])-$l;

  for(; $r < $remaining; $r++) {
    my @bp = map {substr($_,-$r-1,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  return ($l,$r);
}