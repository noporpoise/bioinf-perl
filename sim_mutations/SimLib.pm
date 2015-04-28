package SimLib;

use strict;
use warnings;

use List::Util qw(first min max);
use POSIX;
use Carp;

use GeneticsModule;
use FASTNFile;

# All methods are object methods except these:
use base 'Exporter';
our @EXPORT = qw(normalise_variant print_variant
                 trim_alleles get_flank_shifts get_var_kmers get_key
                 load_genome_mask_files);

sub print_variant
{
  my ($lflank,$rflank,@alleles) = @_;
  print "  LF:$lflank RF:$rflank\n";
  print "  alleles: ".join(',',@alleles)."\n";
}

sub normalise_variant
{
  my ($kmer_size,$lflank,$rflank,@alleles) = @_;

  # 1) Trim alleles
  ($lflank, $rflank, @alleles) = trim_alleles($kmer_size, $lflank, $rflank,
                                              @alleles);

  # 2) Get shift-ability
  my ($shift_left,$shift_right) = get_flank_shifts($lflank, $rflank, @alleles);

  # 2) get kmer keys
  my ($kmer0, $kmer1) = get_var_kmers($kmer_size, $lflank, $rflank,
                                      $shift_left, $shift_right);

  my ($key0, $key1) = map {get_key($_)} ($kmer0, $kmer1);

  # 3) revcmp if needed
  if($key0 gt $key1 || ($key0 eq $key1 && $key0 ne $kmer0 && $key1 eq $kmer1))
  {
    # Reverse complement this variant
    ($lflank, $rflank) = (rev_comp($rflank), rev_comp($lflank));
    @alleles = map {rev_comp($_)} @alleles;
    # next line not needed
    # ($shift_left, $shift_right) = ($shift_right,$shift_left);

    # Shift right if needed
    if($shift_right > 0)
    {
      my $append = substr($rflank, 0, $shift_right);
      substr($rflank, 0, $shift_right) = '';
      $lflank .= $append;

      for(my $i = 0; $i < @alleles; $i++) {
        my $len = length($alleles[$i]);
        if($len > 0) {
          my $dist = $shift_right;
          if($len < $shift_right) { $dist %= $len; }
          substr($alleles[$i], 0, $dist) = '';
          $alleles[$i] .= $append;
        }
      }
    }
  }

  # Remove dupes in alleles
  my %ahsh = ();
  @ahsh{@alleles} = 1;

  @alleles = sort(keys %ahsh);

  return ($lflank, $rflank, @alleles);
}

sub trim_alleles
{
  my ($kmer_size,$lf,$rf,@a) = @_;
  my ($lend,$rend) = get_trim_dist(@a);
  if($lend > 0) {
    @a = map {substr($_, $lend)} @a;
    $lf .= substr($a[0],0,$lend);
    $lf = substr($lf, -$kmer_size);
  }
  if($rend > 0) {
    @a = map {substr($_, 0, -$rend)} @a;
    $rf = substr($a[0],-$rend).$rf;
    $rf = substr($lf,0,$kmer_size);
  }
  return ($lf,$rf,@a);
}

# Get number of matching bases on left and right
sub get_trim_dist
{
  my @a = @_;
  my ($l, $r) = (0, 0);

  my $len = min(map {length($_)} @a);

  for(; $l < $len; $l++) {
    my @bp = map {substr($_,$l,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  my $remaining = $len-$l;

  for(; $r < $remaining; $r++) {
    my @bp = map {substr($_,-$r-1,1)} @a;
    if(defined(first {$_ ne $bp[0]} @bp)) { last; }
  }
  return ($l,$r);
}

sub get_flank_shifts
{
  my ($lf,$rf,@alleles) = @_;
  my ($min_left,$min_right) = (LONG_MAX, LONG_MAX);

  for my $allele (@alleles)
  {
    my $l = length($allele);
    if($l > 0) {
      # left
      my $join = $lf.$allele;
      my $m = 0;
      while($m < length($lf) && substr($join,-$m-1,1) eq substr($join,-$m-$l-1,1)) {
        $m++;
      }
      $min_left = min($min_left, $m);

      # right
      $join = $allele.$rf;
      $m = 0;
      while($m < length($rf) && substr($join,$m,1) eq substr($join,$m+$l,1)) {
        $m++;
      }
      $min_right = min($min_right, $m);
    }
  }
  return ($min_left, $min_right);
}

# Variant key is defined by <key1><key2>
# where <key1> is less than <key2>
# and the keys are the kmer keys of the start of the flanks
# if <key1> == <key2>, then the variant is `forward` when <kmer1>==<key1>
sub get_var_kmers
{
  my ($kmer_size,$lf,$rf,$shift_left,$shift_right) = @_;

  my $joint = $lf.$rf;
  my $kmer0 = substr($joint, length($lf)+$shift_right-$kmer_size, $kmer_size);
  my $kmer1 = substr($joint, -length($rf)-$shift_left, $kmer_size);

  return ($kmer0, $kmer1);
}

sub get_key
{
  my ($kmer) = @_;
  my $rev_comp = rev_comp($kmer);
  return ($rev_comp lt $kmer ? $rev_comp : $kmer);
}

# All masks and chromosomes should be the same length
sub load_genome_mask_files
{
  my @files = @_;
  my @genomes = ();
  my @masks = ();
  my ($chrname,$alnlen,$title,$seq);
  my $fastn;

  while(@files > 0)
  {
    my $genome_file = shift(@files);
    my $mask_file = shift(@files);
    for my $file ($genome_file, $mask_file) {
      $fastn = open_fastn_file($file);
      ($title,$seq) = $fastn->read_next();
      close_fastn_file($fastn);
      if(!defined($title)) { die("Empty file: $file\n"); }
      if(!defined($chrname)) { $chrname = $title; $alnlen = length($seq); }
      elsif(length($seq) != $alnlen) {
        die("Genomes diff lengths [file: $file $alnlen vs ".length($seq)."]");
      }
      if($file eq $genome_file) { push(@genomes, uc($seq)); }
      else { push(@masks, $seq); }
    }
  }

  # Remove fasta/fastq 'comments' (only take chars upto first whitespace)
  $chrname =~ s/\s.*$//g;
  return (\@genomes, \@masks, $chrname, $alnlen);
}

1;
