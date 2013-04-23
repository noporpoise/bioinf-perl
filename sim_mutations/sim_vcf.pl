#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first min);
use File::Path qw(make_path);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use FASTNFile;
use GeneticsModule;
use UsefulModule;
use SimLib;

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
# print STDERR "$mask\n";

#
# Generate truth VCF
#
print "##".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT))."\n";
my $prev_flank_len = -1;

for(my $i = 0; $mask =~ /([^\.]+(?:\.+[^\.]+)*?)(\.{$kmer_size,1000})/g; $i++)
{
  # $-[0] is start pos of allele
  my $start = $-[0];
  my $var_len = length($1);
  my $mutations = $1;

  my $lflank_len = $prev_flank_len;
  my $rflank_len = length($2);

  if($lflank_len == -1)
  {
    $lflank_len = 0;
    while($lflank_len < $start && $lflank_len < 1000 &&
          substr($mask, $start-$lflank_len-1, 1) eq '.') { $lflank_len++; }
  }
  $prev_flank_len = $rflank_len;

  # print STDERR "Match '$1' start:$start len:$var_len mut:$mutations\n";

  # Get alleles
  my @alleles = map {substr($genomes[$_], $start, $var_len)} @sampleids;
  for(my $i = 0; $i < @alleles; $i++) { $alleles[$i] =~ s/-//g; }

  # Remove duplicates
  my %ah = ();
  @ah{@alleles} = 1;
  @alleles = keys(%ah);

  if(@alleles > 1 && $lflank_len >= $kmer_size)
  {
    my $lflank = substr($genomes[0], $start-$lflank_len, $lflank_len);
    my $rflank = substr($genomes[0], $start+$var_len, $rflank_len);

    ($lflank, $rflank, @alleles) = normalise_variant($kmer_size, $lflank, $rflank,
                                                     @alleles);

    # If not a SNP add a padding base from the left flank
    if(defined(first {length($_) != 1} @alleles)) {
      map {$alleles[$_] = substr($lflank,-1).$alleles[$_]} 0..$#alleles;
    }

    my $info = "LF=".substr($lflank,-$kmer_size).";" .
               "RF=".substr($rflank,0,$kmer_size);

    # Get annotations from the mask
    if($mutations =~ /S/) { $info .= ";SNP"; }
    if($mutations =~ /I/i) { $info .= ";INS"; }
    if($mutations =~ /D/i) { $info .= ";DEL"; }
    if($mutations =~ /V/) { $info .= ";INV"; }

    print join("\t", ".", "0", "var$i", "N", join(',', @alleles), ".", "PASS",
               $info, ".") . "\n";
  }
}