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

#
# Load and merge mask files
#
my $num_of_samples = int(@ARGV / 2);

# All masks and chromosomes should be the same length

my @genomes = ();
my @masks = ();
my @sampleids = (0..($num_of_samples-1));

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

print STDERR "Genomes and masks loaded\n";

#
# Merge masks
#
my $mask = "." x length($masks[0]);

for(my $i = 0; $i < $num_of_samples; $i++)
{
  while($masks[$i] =~ /([^\.]+)/g)
  {
    my $start = $-[0];
    my $end = $+[0];
    substr($mask, $start, $end-$start) = $1;
  }
}

# print ">mask\n$mask\n";

print STDERR "Masks overlaid\n";

#
# Generate truth VCF
#
print "##".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT))."\n";

for(my $i = 0; $mask =~ /([^\.]+(?:\.+[^\.]+)*?)\.{$kmer_size}/g; $i++)
{
  # $-[0] is start pos of allele
  my $start = $-[0];
  my $len = length($1);

  # print "Match '$1' start:$start len:$len\n";

  # Get alleles
  my @alleles = map {substr($genomes[$_], $start, $len)} @sampleids;
  map {$_ =~ s/\-//g} @alleles;

  # Remove duplicates
  my %ah = ();
  @ah{@alleles} = 1;
  @alleles = keys(%ah);

  if(@alleles > 1 && $start >= $kmer_size)
  {
    # Trim matching bp off ends
    my ($lend,$rend) = trim_alleles(@alleles);

    if($lend > 0 || $rend > 0) {
      $start += $lend;
      $len -= $lend + $rend;
      @alleles = map {substr($_, $lend, -$rend)} @alleles;
    }

    my $lflank = substr($genomes[0], $start-$kmer_size, $kmer_size);
    my $rflank = substr($genomes[0], $start+$len, $kmer_size);

    if(defined(first {length($_) != 1} @alleles))
    {
      # Not a SNP
      map {$alleles[$_] = substr($lflank,-1).$alleles[$_]} 0..$#alleles;
    }

    my $info = "LF=$lflank;RF=$rflank";

    # Get annotations from the mask
    my $m = substr($mask, $start, $len);
    if($m =~ /S/) { $info .= ";SNP"; }
    if($m =~ /I/i) { $info .= ";INS"; }
    if($m =~ /D/i) { $info .= ";DEL"; }
    if($m =~ /V/) { $info .= ";INV"; }

    print join("\t", ".", "0", "var$i", "N", join(',', @alleles), ".", "PASS",
               $info, ".") . "\n";
  }
}

# Get number of matching bases on left and right
sub trim_alleles
{
  my @a = @_;
  my ($l, $r) = (0, 0);
  # print "Trim:".join(',', @a)."\n";

  for(; $l < length($_[0]); $l++) {
    my @bp = map {substr($_,$l,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  my $remaining = length($_[0])-$l;

  for(; $r < $remaining; $r++) {
    my @bp = map {substr($_,-$r-1,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  return ($l,$r);
}