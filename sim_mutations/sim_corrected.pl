#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max reduce);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use GeneticsModule;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./sim_corrected.pl [--print-mismatches] <file.fa> <ref.fa ..>
  Prints stats on read match/mismatch rate.
\n";
  exit(-1);
}

my $print_mismatches = 0;

while(@ARGV > 0 && $ARGV[0] =~ /^--/) {
  if($ARGV[0] =~ /^--print-mismatche?s?$/) {
    shift;
    $print_mismatches = 1;
  } else {
    print_usage("Unknown option: $ARGV[0]");
  }
}

if(@ARGV < 2) { print_usage(); }

my $reads_path = shift(@ARGV);
my @ref_paths = @ARGV;

# Open reads file
my $fastn = open_fastn_file($reads_path);

# 1) Load ref genome

print STDERR "Loading ref...\n";
my ($genome) = read_all_from_files(@ref_paths);

my $genome_size = 0;

while(my ($chrom,$chr_seq) = each(%$genome)) {
  print STDERR "Chrom: $chrom\n";
  # $genome->{$chrom} = uc($chr_seq);
  $genome_size += length($chr_seq);
}

print STDERR "Genome size: ".num2str($genome_size)."\n";
print STDERR "Loading reads...\n";

# Stats
my ($base_count, $base_match, $base_mismatch) = (0,0,0);
my ($uc_count,   $uc_match,   $uc_mismatch)   = (0,0,0);
my ($lc_count,   $lc_match,   $lc_mismatch)   = (0,0,0);
my $num_Ns = 0;
my $num_reads++;

my ($orig_len_sum, $orig_len_diff) = (0,0);
my ($changed_total, $changed_match, $changed_mismatch) = (0,0,0,0);

# 2) Parse reads
my ($start,$truth);
my ($name,$seq);
while((($name,$seq) = $fastn->read_next()) && defined($name))
{
  # name is >r2692:genome0:47174:47693/1 orig=ACAATCACAGAGATTGTAAGT
  if($name =~ /^.*?:(.*?):(\d+):(\d+)(\/[12])?(?:.*?orig=(\w+))?/) {
    my ($chrom,$startr1,$startr2,$rnum,$orig_seq) = ($1,$2,$3,$4,$5);
    if(!defined($rnum)) { $rnum = '/1'; }
    if($rnum eq '/1') { $start = $startr1; }
    else              { $start = $startr2; }
    if(!defined($genome->{$chrom})) { die("Cannot find chrom: $chrom"); }
    $truth = substr($genome->{$chrom}, $start, length($seq));

    # Compare $seq to $truth
    # $orig_seq may not be defined, as it is not always in the output
    compare_reads($name, $seq, $truth, $orig_seq);
  }
  else { die("Bad read: $name"); }
  $num_reads++;
}

# Calc total stats
$base_mismatch = $uc_mismatch + $lc_mismatch;
$base_match = $base_count - $base_mismatch;

# Calc uppercase stats
$uc_match = $uc_count - $uc_mismatch;

# Calc lowercase stats
$lc_match = $lc_count - $lc_mismatch;

# Print stats
print STDERR "Uppercase:\n";
print STDERR "     total: ".pretty_fraction($uc_count,   $base_count)."\n";
print STDERR "     match: ".pretty_fraction($uc_match,   $uc_count)  ."\n";
print STDERR "  mismatch: ".pretty_fraction($uc_mismatch,$uc_count)  ."\n";
print STDERR "Lowercase:\n";
print STDERR "     total: ".pretty_fraction($lc_count,   $base_count)."\n";
print STDERR "     match: ".pretty_fraction($lc_match,   $lc_count)  ."\n";
print STDERR "  mismatch: ".pretty_fraction($lc_mismatch,$lc_count)  ."\n";
print STDERR "All:\n";
print STDERR "     match: ".pretty_fraction($base_match,   $base_count)."\n";
print STDERR "  mismatch: ".pretty_fraction($base_mismatch,$base_count)."\n";
print STDERR "   N bases: ".pretty_fraction($num_Ns,       $base_count)."\n";
print STDERR "     reads: ".num2str($num_reads)."\n";

if($orig_len_sum > 0) {
  print STDERR "Bases Changed:\n";
  print STDERR "     change info: ".pretty_fraction($orig_len_sum,     $base_count)   ."\n";
  print STDERR "   bases changed: ".pretty_fraction($changed_total,    $orig_len_sum) ."\n";
  print STDERR "    change match: ".pretty_fraction($changed_match,    $changed_total)."\n";
  print STDERR " change mismatch: ".pretty_fraction($changed_mismatch, $changed_total)."\n";
  print STDERR "   read len diff: ".pretty_fraction($orig_len_diff,    $base_count)   .
               " (not included in mismatch rate)\n";
}

print STDERR "Coverage:\n";
print_covgerage('ACGTN coverage', $base_count,         $genome_size);
print_covgerage('ACGT  coverage', $base_count-$num_Ns, $genome_size);

close_fastn_file($fastn);

sub compare_reads
{
  my ($title,$read,$ref,$orig_seq) = @_;
  my $len = length($read);
  my $read_Ns = 0;

  # If we have the orignal sequence, gather simple correction stats
  if(defined($orig_seq)) {
    my $orig_len = length($orig_seq);
    $orig_len_sum += $orig_len;
    my $minlen = min($orig_len, $len);
    my $maxlen = max($orig_len, $len);
    for(my $i = 0; $i < $minlen; $i++) {
      my ($o,$r,$t) = map {uc(substr($_,$i,1))} ($orig_seq,$read,$truth);
      if($o ne $r) {
        $changed_total++;
        if($r eq $t) { $changed_match++;    }
        else         { $changed_mismatch++; }
      }
    }
    $orig_len_diff += $maxlen - $minlen;
  }

  my $prev_mismatches = $uc_mismatch + $lc_mismatch;

  for(my $i = 0; $i < $len; $i++) {
    my ($a,$b) = map {substr($_,$i,1)} ($read,$ref);
    if(uc($a) eq 'N') { $read_Ns++; }
    else {
      my $isupper = (uc($a) eq $a);
      if(uc($a) ne uc($b)) {
        $uc_mismatch +=  $isupper;
        $lc_mismatch += !$isupper;
      }
      $uc_count +=  $isupper;
      $lc_count += !$isupper;
    }
  }

  if($print_mismatches && $uc_mismatch + $lc_mismatch > $prev_mismatches) {
    print ">$title\n$read\n$truth\n";
  }

  $base_count += $len;
  $num_Ns += $read_Ns;
}

sub print_covgerage
{
  my ($title,$nom,$denom) = @_;
  print STDERR "  $title: ".sprintf("%.1fX",$nom / $denom)." " .
               "(" . num2str($nom) . "/" . num2str($denom) . ")\n";
}
