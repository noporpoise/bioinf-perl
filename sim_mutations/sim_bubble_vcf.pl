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
"Usage: ./sim_bubble_vcf.pl <kmer_size> <in1.fa> <in1.mask> <in2.fa> <in2.mask>
  Print bubble vcf\n";

  exit(-1);
}

if(@ARGV != 5) { print_usage(); }

my $kmer_size = shift;

if($kmer_size !~ /^\d+$/) { print_usage(); }
if(scalar(@ARGV) % 2 != 0) { print_usage(); }

#
# Load and merge mask files
#
my ($chrname,$len,$genarr,$mskarr) = load_genome_mask_files(@ARGV);
my @seq = @$genarr;
my @masks = @$mskarr;

print STDERR "".@seq." Genome and mask pairs loaded\n";

#
# Generate truth VCF
#
print "#".join("\t", qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT))."\n";

# Find position of first mismatch
my $clamp_start = next_clamp(0);

if($clamp_start == -1) {
  print STDERR "Not a single stretch of identity of length >=$kmer_size\n";
  exit;
}

my ($ls,$le,$rs,$re);

for(my $n = 0; (($le,$rs) = next_bubble($clamp_start)) && defined($le); $n++)
{
  my ($ls, $re) = (extend_left_clamp($le), extend_right_clamp($rs));
  my $lf = substr($seq[0],$ls,$le-$ls+1);
  my $rf = substr($seq[0],$rs,$re-$rs+1);
  $lf =~ s/-//g; $rf =~ s/-//g;

  my @branches = map{'N'.substr($seq[$_], $le+1, $rs-$le-1)} (0,1);
  map {$branches[$_] =~ s/-//g} (0,1);

  my @m = map{substr($masks[$_], $le+1, $rs-$le-1)} (0,1);
  my ($snp,$ins,$del,$inv) = mutant_breakdown(@m);

  my $info = "LF=$lf;RF=$rf;SNP=$snp;INS=$ins;DEL=$del;INV=$inv";

  print join("\t", $chrname, "0", "truth$n", "N", join(',', @branches), ".",
             "PASS", $info, ".") . "\n";

  $clamp_start = $rs;
}

# pass two mask strings
# returns (snps,ins,del,inv) counts
sub mutant_breakdown
{
  my @m = @_;
  my ($i, $mlen) = (0, length($m[0]));
  my %vars = ('S' => 0, 'I' => 0, 'D' => 0, 'V' => 0);
  while($i < $mlen) {
    my @v = map{uc(substr($m[$_],$i,1))} (0,1);
    if($v[0] ne $v[1]) {
      my $j = $v[0] eq '.' ? 1 : 0;
      my $var = $v[$j];
      $vars{$var}++;
      do { $i++; } while($i < $mlen && substr($m[$j],$i,1) eq $var);
    } else { $i++; }
  }
  return @vars{qw(S I D V)};
}

# returns position of matching base before and base after variants
sub next_bubble
{
  my ($clend, $clstart) = (0, @_);

  # Extend left-hand match until first mismatch
  for($clend = $clstart+$kmer_size;
      $clend < $len && substr($seq[0],$clend,1) eq substr($seq[1],$clend,1);
      $clend++) {}

  if($clend == $len) { return undef; }

  # Find next matching run of $kmer_size bp
  my $pos = next_clamp($clend+1);

  return ($pos == -1) ? undef : ($clend-1, $pos);
}

# A clamp is just a run of $kmer_size bases of matching sequence
sub next_clamp
{
  my ($c, $match, $run, $pos) = (0, 0, 0, @_);
  while($match < $kmer_size && $pos < $len)
  {
    $c = substr($seq[0],$pos,1);
    if($c eq substr($seq[1],$pos,1)) {
      if($c ne '-') { $match++; }
      $run++;
    }
    else { $match = $run = 0; }
    $pos++;
    # print "$c: pos: $pos match: $match run: $run\n";
  }
  return ($match < $kmer_size ? -1 : $pos-$run);
}

sub extend_right_clamp
{
  my ($match, $p) = (1, @_);
  while($match < $kmer_size) {
    if(substr($seq[0],$p,1) ne '-') { $match++; }
    $p++;
  }
  return $p;
}

sub extend_left_clamp
{
  my ($match, $p) = (1, @_);
  while($match < $kmer_size) {
    if(substr($seq[0],$p-1,1) ne '-') { $match++; }
    $p--;
  }
  return $p;
}
