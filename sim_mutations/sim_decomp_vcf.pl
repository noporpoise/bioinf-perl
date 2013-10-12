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
my $ref_path = shift;

if(@ARGV % 2 != 0) { print_usage("Expected odd number of args"); }

#
# Load and merge mask files
#
# All masks and chromosomes should be the same length

my @genomes = ();
my @masks = ();
my ($chrname,$ref);
my ($title,$seq);
my $fastn;

$fastn = open_fastn_file($ref_path);
($chrname,$ref) = $fastn->read_next();
close_fastn_file($fastn);
if(!defined($chrname)) { die("Empty file: $ref_path\n"); }
$ref = uc($ref);

# Remove fasta/fastq 'comments' (only take chars upto first whitespace)
$chrname =~ s/\s.*$//g;

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
print "##reference=file://$ref_path\n";
print "##FILTER=<ID=PASS,Description=\"All filters passed\"\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##contig=<ID=un,length=1000000,assembly=None>\n";
print "##contig=<ID=$chrname,length=$len>\n";
print "#".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE))."\n";

my ($start, $end, $sample, $i) = (0,0,-1);
my @alleles;
my %hsh;
for(my $var = 0; ; $var++)
{
  ($start,$end,$sample) = get_var($start, $sample+1);
  if(!defined($start)) { last; }
  if($start == 0) { next; }

  # Get alleles, remove deletions ('-')
  @alleles = map {substr($_, $start, $end-$start)} @genomes;
  map {$alleles[$_] =~ s/\-//g} 0..$#alleles;

  my $r = substr($ref, $start, $end-$start);

  # Remove duplicates
  %hsh = ();
  @hsh{@alleles} = 1;
  delete($hsh{$r}); # remove ref allele from alts
  @alleles = keys(%hsh);

  # Add padding base
  my $pos = $start;
  if(defined(first {length($_) != 1} @alleles)) {
    $pos--;
    my $c = substr($ref, $pos, 1);
    ($r, @alleles) = map {$c.$_} ($r, @alleles);
  }

  my $alt = join(',', @alleles);
  print join("\t", $chrname, $pos+1, "truth$var", $r, $alt, '.', "PASS",
             ".", "GT", "0/1")."\n";
}

sub get_var
{
  # In mask files, variants start with lower case letters
  # s = SNPs; i = insert; d = deletion; v = inversion
  # see sim_mutation.pl
  my ($pos,$sample) = @_;
  my ($i,$c,$end);
  if($sample == @masks) { $sample = 0; $pos++; }
  for(; $pos < $len; $pos++) {
    for($i = $sample; $i < @masks; $i++) {
      $c = substr($masks[$i],$pos,1);
      if($c ne '.' && $c eq lc($c)) {
        $c = uc($c);
        for($end = $pos+1; $end < $len && substr($masks[$sample],$end,1) eq $c; $end++) {}
        return ($pos,$end,$sample);
      }
    }
    $sample = 0;
  }
  return ();
}
