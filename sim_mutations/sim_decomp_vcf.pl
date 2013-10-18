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
"Usage: ./sim_decomp_vcf.pl [<sample1.fa> <sample1.mask> ..]\n";

  exit(-1);
}

if(@ARGV == 0) { print_usage(); }
if(@ARGV % 2 != 0) { print_usage("Expected odd number of args"); }

my $ref_path = $ARGV[0];

#
# Load and merge mask files
#
my ($chrname,$len,$genarr,$mskarr) = load_genome_mask_files(@ARGV);
my @genomes = @$genarr;
my @masks = @$mskarr;

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

my ($start, $end, $refpos) = (0,0,0);
my $c;
my @alleles;
my %hsh;
for(my $var = 0; ; $var++)
{
  ($start,$end,$refpos) = get_var($end, $refpos);
  if(!defined($start)) { last; }
  if($start == 0) { next; }

  # Get alleles, remove deletions ('-')
  @alleles = map {substr($_, $start, $end-$start)} @genomes;
  map {$alleles[$_] =~ s/\-//g} 0..$#alleles;

  my @m = map {substr($masks[$_], $start, $end-$start)} 0..$#masks;

  my $r = $alleles[0];

  # Remove duplicates
  %hsh = ();
  @hsh{@alleles} = 1;
  delete($hsh{$r}); # remove ref allele from alts
  @alleles = keys(%hsh);

  if(@alleles == 0) { next; }

  # Add padding base
  my $varpos = $refpos;
  if(defined(first {length($_) != 1} ($r,@alleles))) {
    $varpos--;
    my $pos = $start-1;
    while($pos >= 0 && ($c = substr($genomes[0], $pos, 1)) eq '-') { $pos--; }
    if($pos < 0) { $c = 'N'; }
    ($r, @alleles) = map {$c.$_} ($r, @alleles);
  }

  my $alt = join(',', @alleles);
  my $info = "."; # "L=$start:$end;D=".join(',', @m);
  print join("\t", $chrname, $varpos+1, "truth$var", $r, $alt, '.', "PASS",
             $info, "GT", "0/1")."\n";

  # Move refpos forward
  for(; $start < $end; $start++) {
    if(substr($genomes[0],$start,1) ne '-') { $refpos++; }
  }
}

sub get_var
{
  # In mask files, variants start with lower case letters
  # s = SNPs; i = insert; d = deletion; v = inversion
  # Assume no overlapping variants; see sim_mutation.pl
  my ($pos,$refpos) = @_;
  my ($m,$c,$end);
  for(; $pos < $len; $pos++) {
    for($m = 0; $m < @masks; $m++) {
      $c = substr($masks[$m],$pos,1);
      if($c ne '.' && $c eq lc($c)) {
        $c = uc($c);
        for($end = $pos+1; $end < $len && substr($masks[$m],$end,1) eq $c; $end++) {}
        return ($pos,$end,$refpos);
      }
    }
    if(substr($genomes[0],$pos,1) ne '-') { $refpos++; }
  }
  return ();
}
