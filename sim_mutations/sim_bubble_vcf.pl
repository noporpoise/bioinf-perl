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
"Usage: ./sim_bubble_vcf.pl <kmer_size> [<sample.fa> <sample.mask> ..]\n";

  exit(-1);
}

if(@ARGV < 3) { print_usage(); }

my $kmer_size = shift;

if($kmer_size !~ /^\d+$/) { print_usage(); }
if(scalar(@ARGV) % 2 != 0) { print_usage(); }

#
# Load and merge mask files
#
my ($chrname,$len,$genarr,$mskarr) = load_genome_mask_files(@ARGV);
my @genomes = @$genarr;
my @masks = @$mskarr;

print STDERR "".@genomes." Genome and mask pairs loaded\n";

#
# Merge masks
#
my $mask = "." x length($masks[0]);

for(my $i = 0; $i < @genomes; $i++)
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
print "#".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT))."\n";
my $prev_flank_len = -1;

# speed up
if($mask =~ /^\.*$/) { print STDERR "No variants\n"; exit; }

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
  my @alleles = map {substr($genomes[$_], $start, $var_len)} 0..$#genomes;
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
    my @local_masks = map {substr($masks[$_],$start,$var_len)} 0..$#genomes;

    # The start of each variant is lower case in only one individual
    my ($snp,$ins,$del,$inv) = (0,0,0,0);
    for my $local_mask (@local_masks)
    {
      while($local_mask =~ /s/g) { $snp++; }
      while($local_mask =~ /i/g) { $ins++; }
      while($local_mask =~ /d/g) { $del++; }
      while($local_mask =~ /v/g) { $inv++; }
    }

    $info .= ";SNP=$snp;INS=$ins;DEL=$del;INV=$inv";

    print join("\t", $chrname, "0", "truth$i", "N", join(',', @alleles), ".",
               "PASS", $info, ".") . "\n";
  }
}
