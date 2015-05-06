#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(first min);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

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
"Usage: ./sim_vcf.pl <ref.fa> [<sample1.fa> <sample1.mask> ..]\n";

  exit(-1);
}

my $command = "$0 @ARGV";
# my $user_chrname;

# while(@ARGV && $ARGV[0] =~ /^-/) {
#   my $cmd = shift;
#   if($cmd eq "-n") { $user_chrname = shift; }
#   else { print_usage("Unknown command: $cmd"); }
# }

if(@ARGV == 0) { print_usage(); }
if(@ARGV % 2 != 1) { print_usage("Expected odd number of args"); }

my $ref = shift;
my @files = @ARGV;

#
# Load and merge mask files
#
my $fastn = open_fastn_file($ref);
my ($ref_chrom, $ref_seq) = $fastn->read_next();
if(!defined($ref_seq)) { die("Couldn't read a sequence from ref: $ref"); }
close_fastn_file($fastn);

$ref_seq = uc($ref_seq);
# Remove comments after title (everything after first whitespace)
($ref_chrom) = ($ref_chrom =~ /^(\S*)/);
print STDERR "ref: '$ref_chrom'\n";

my ($genarr,$mskarr,undef,$chrlen) = load_genome_mask_files(@files);
my @genomes = @$genarr;
my @masks = @$mskarr;

if($chrlen != length($ref_seq)) { die("Ref seq len mismatch"); }

# if(defined($user_chrname)) { $chrname = $user_chrname; }

print STDERR "".@genomes." Genome and mask pairs loaded\n";

# map {print "$_: ".$genomes[$_]."\n"; } 0..$#genomes;

#
# Generate truth VCF
#
print "##fileformat=VCFv4.1\n";
print "##fileDate=20130930\n";
print "##reference=unknown\n";
print "##cmd=$ref\n";
print "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print "##contig=<ID=un,length=1000000,assembly=None>\n";
print "##contig=<ID=$ref_chrom,length=$chrlen>\n";
print "#".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE))."\n";

my ($start, $end) = (0,0,0);
my $c;
my @alleles;
my %hsh;

for(my $var = 0; ; $var++)
{
  ($start,$end) = next_var($end,$chrlen);
  if(!defined($start)) { last; }
  # if($start == 0) { next; } # Don't report events at the very start of the chrom

  # Get alleles, remove deletions ('-')
  @alleles = map {substr($_, $start, $end-$start)} @genomes;
  map {$alleles[$_] =~ s/\-//g} 0..$#alleles;

  # my @m = map {substr($masks[$_], $start, $end-$start)} 0..$#masks;

  my $r = substr($ref_seq, $start, $end-$start);

  # Remove duplicates
  %hsh = ();
  @hsh{@alleles} = 1;
  delete($hsh{$r}); # remove ref allele from alts
  @alleles = keys(%hsh);

  if(@alleles == 0) { next; }

  # Add padding base
  my $varpos = $start;
  if(defined(first {length($_) != 1} ($r,@alleles))) {
    $varpos--;
    my $pos = $start-1;
    while($pos >= 0 && ($c = substr($genomes[0], $pos, 1)) eq '-') { $pos--; }
    if($pos < 0) { $c = 'N'; }
    ($r, @alleles) = map {$c.$_} ($r, @alleles);
  }

  my $alt = join(',', @alleles);
  my $info = ".";
  # my $info = "L=$start:$end;D=".join(',', @m);

  print join("\t", $ref_chrom, $varpos+1, "truth$var", $r, $alt, '.', "PASS",
             $info, "GT", "./.")."\n";
}

sub next_var
{
  # In mask files, variants start with lower case letters
  # s = SNPs; i = insert; d = deletion; v = inversion
  # Assume no overlapping variants; see sim_mutation.pl
  my ($pos,$len) = @_;
  my ($m,$n,$c,$end);

  for(; $pos < $len; $pos++)
  {
    for($m = 0; $m < @masks; $m++)
    {
      $c = substr($masks[$m],$pos,1);

      if($c ne '.' && $c eq lc($c) &&
         defined(first {substr($masks[$_],$pos,1) ne $c} 0..$#masks))
      {
        $c = uc($c);
        for($end = $pos+1; $end < $len && substr($masks[$m],$end,1) eq $c; $end++) {}
        return ($pos,$end);
      }
    }
  }
  return ();
}
