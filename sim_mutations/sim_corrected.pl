#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(reduce);

# Use current directory to find modules
# use FindBin;
# use lib ../lib;

use FASTNFile;
use GeneticsModule;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./sim_corrected.pl <file.fa> <ref.fa ..>
  Prints stats on read match/mismatch rate.
\n";
  exit(-1);
}

if(@ARGV < 2) { print_usage(); }

my $reads_path = shift(@ARGV);
my @ref_paths = @ARGV;

# Open reads file
my $fastn = open_fastn_file($reads_path);

# 1) Load ref genome

print "Loading ref...\n";
my ($genome) = read_all_from_files(@ref_paths);

my $genome_size = 0;

while(my ($chrom,$chr_seq) = each(%$genome)) {
  print "Chrom: $chrom\n";
  $genome_size += length($chr_seq);
}

print "Genome size: ".num2str($genome_size)."\n";
print "Loading reads...\n";

# Stats
my ($base_count, $base_match, $base_mismatch) = (0,0,0);
my ($uc_count,   $uc_match,   $uc_mismatch)   = (0,0,0);
my ($lc_count,   $lc_match,   $lc_mismatch)   = (0,0,0);
my $num_Ns = 0;
my $num_reads++;

# 2) Parse reads
my ($start,$truth);
my ($name,$seq);
while((($name,$seq) = $fastn->read_next()) && defined($name))
{
  # name is >r2692:genome0:47174:47693/1
  if($name =~ /^.*?:(.*?):(\d+):(\d+)(\/[12])?/) {
    my ($chrom,$startr1,$startr2,$rnum) = ($1,$2,$3,$4);
    if(!defined($rnum)) { $rnum = '/1'; }
    if($rnum eq '/1') { $start = $startr1; }
    else              { $start = $startr2; }
    if(!defined($genome->{$chrom})) { die("Cannot find chrom: $chrom"); }
    $truth = substr($genome->{$chrom}, $start, length($seq));

    # Compare $seq to $truth
    compare_reads($seq, $truth);
  }
  else { die("Bad read: $name"); }
  $num_reads++;
}

# Calc total stats
$base_match = $base_count - $base_mismatch;

# Calc uppercase stats
$uc_match = $uc_count - $uc_mismatch;

# Calculate lowercase stats
$lc_count = $base_count - $uc_count;
$lc_mismatch = $base_mismatch - $uc_mismatch;
$lc_match = $base_match - $uc_match;

# DEV: add number of reads changed

# Print stats
print "Uppercase:\n";
print "     total: ".pretty_fraction($uc_count,   $base_count)."\n";
print "     match: ".pretty_fraction($uc_match,   $uc_count)."\n";
print "  mismatch: ".pretty_fraction($uc_mismatch,$uc_count)."\n";
print "Lowercase:\n";
print "     total: ".pretty_fraction($lc_count,   $base_count)."\n";
print "     match: ".pretty_fraction($lc_match,   $lc_count)."\n";
print "  mismatch: ".pretty_fraction($lc_mismatch,$lc_count)."\n";
print "All:\n";
print "     match: ".pretty_fraction($base_match,   $base_count)."\n";
print "  mismatch: ".pretty_fraction($base_mismatch,$base_count)."\n";
print "   N bases: ".pretty_fraction($num_Ns,       $base_count)."\n";
print "     reads: ".num2str($num_reads)."\n";
print_covgerage('ACGTN coverage', $base_count,         $genome_size);
print_covgerage('ACGT  coverage', $base_count-$num_Ns, $genome_size);

close_fastn_file($fastn);

sub compare_reads
{
  my ($read,$ref) = @_;
  my $len = length($read);
  my $read_Ns = 0;

  for(my $i = 0; $i < $len; $i++) {
    my ($a,$b) = map {substr($_,$i,1)} ($read,$ref);
    if(uc($a) eq 'N') { $read_Ns++; }
    else {
      my $isupper = (uc($a) eq $a);
      if(uc($a) ne uc($b)) {
        if($isupper) { $uc_mismatch++; }
        $base_mismatch++;
      }
      $uc_count += $isupper;
    }
  }

  $base_count += $len;
  $num_Ns += $read_Ns;
}

sub print_covgerage
{
  my ($title,$nom,$denom) = @_;
  print "  $title: ".sprintf("%.1fX",$nom / $denom)." " .
        "(" . num2str($nom) . "/" . num2str($denom) . ")\n";
}
