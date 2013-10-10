#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first min);

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
"Usage: ./sim_decomp_vcf.pl <ref.fa> [<sample1.fa> <sample1.mask> ..]\n";

  exit(-1);
}

if(@ARGV < 3) { print_usage(); }
my $ref_file = shift;

if(@ARGV % 2 == 0) { print_usage("Expected odd number of args"); }

#
# Load and merge mask files
#
# All masks and chromosomes should be the same length

my @genomes = ();
my @masks = ();
my ($chrname,$ref);
my ($title,$seq);
my $fastn;

$fastn = open_fastn_file($ref_file);
($chrname,$ref) = $fastn->read_next();
close_fastn_file($fastn);
if(!defined($chrname)) { die("Empty file: $ref_file\n"); }
$ref = uc($ref);

my $len = length($ref);

while(@ARGV > 0)
{
  my $genome_file = shift;
  my $mask_file = shift;
  for my $file ($genome_file, $mask_file) {
    $fastn = open_fastn_file($file);
    ($title,$seq) = $fastn->read_next();
    close_fastn_file($fastn);
    if(!defined($title)) { die("Empty file: $file\n"); }
    if(length($seq) != $len)
    { die("Genomes diff lengths [file: $file $len vs ".length($seq)."]"); }
    if($file eq $genome_file) { push(@genomes, uc($seq)); }
    else { push(@masks, $seq); }
  }
}

print STDERR "".@genomes." Genome and mask pairs loaded\n";

# map {print "$_: ".$genomes[$_]."\n"; } 0..$#genomes;

#
# Generate truth VCF
#
print "##fileformat=VCFv4.1\n";
print "##fileDate=20130930\n";
print "##reference=file://$ref_file\n";
print "##FILTER=<ID=PASS,Description=\"All filters passed\"\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##contig=<ID=un,length=1000000,assembly=None>\n";
print "##contig=<ID=$chrname,length=$len>\n";
print "#".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE))."\n";

my ($start, $end, $i) = (0,0);
my @alleles;
my %hsh;
for(my $var = 0; ($start = get_var_start($end)) != -1; $var++)
{
  $end = get_var_end($start);
  if($start == 0) { next; }

  # Get alleles, remove deletions ('-')
  @alleles = map {substr($_, $start, $end-$start)} @genomes;
  map {$alleles[$_] =~ s/\-//g} 0..$#alleles;

  # Remove duplicates, keep reference allele at the front
  my $r = $alleles[0];
  %hsh = ();
  @hsh{@alleles} = 1;
  delete($hsh{$r});
  @alleles = keys(%hsh);
  unshift(@alleles, $r);

  # Add padding base
  if(defined(first {length($_) != 1} @alleles)) {
    $start--;
    my $c = substr($genomes[0], $start, 1);
    @alleles = map {$c.$_} @alleles;
  }

  my $ref = shift(@alleles);
  my $alt = join(',', @alleles);
  print join("\t", $chrname, $start+1, "truth$var", $ref, $alt, '.', "PASS",
             ".", "GT", "0/1")."\n";
}

sub get_var_start
{
  my ($pos) = @_;
  my ($i,$c);
  for(; $pos < $len; $pos++)
  {
    $c = substr($genomes[0],$pos,1);
    if($c eq '-') { return $pos; }
    for($i = 1; $i < @genomes && substr($genomes[$i],$pos,1) eq $c; $i++) {}
    if($i < @genomes) { return $pos; }
  }
  return -1;
}

sub get_var_end
{
  my ($pos) = @_;
  my ($i,$c);
  for($pos++; $pos < $len; $pos++)
  {
    if(($c = substr($genomes[0],$pos,1)) ne '-')
    {
      for($i = 1; $i < @genomes && substr($genomes[$i],$pos,1) eq $c; $i++) {}
      if($i == @genomes) { return $pos; }
    }
  }
  return $len;
}
