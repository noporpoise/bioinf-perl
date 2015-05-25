#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max min first sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use GeneticsModule;
use UsefulModule;

use constant { FORWARD => 0, REVERSE => 1 };

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  
  print STDERR "" .
"Usage: ./sim_breakpoints.pl <N> <ref.fa> <out.fa> <out.txt>\n" .
"  Simulate <N> breakpoints\n";

  exit(-1);
}

if(@ARGV != 4) { print_usage(); }

my ($N, $ref_path, $out_fa_path, $out_txt_path) = @ARGV;

open(FA, ">$out_fa_path") or die("Cannot open $out_fa_path");
open(TXT, ">$out_txt_path") or die("Cannot open $out_txt_path");

# Load sequence
my $fastn = open_fastn_file($ref_path);
my ($title,$seq) = $fastn->read_next();
if(!defined($title)) { die("Woops, can't read sequence: $ref_path"); }
close_fastn_file($fastn);

$seq = uc($seq);
($title) = ($title =~ /^(\S*)/);
print STDERR "ref: '$title'\n";

my $ref_len = length($seq);

# Generate N breakpoints
# each point is the start of a new block
my @points = sort {$a <=> $b} map {int(rand($ref_len-1))} 1..$N;

my @blocks = ();
my $start = 0;
my $min_block = 40;

# Remove blocks too short
print STDERR "Removing blocks shorter than $min_block\n";
my @points_tmp = ();
for(my $i = 0; $i < @points; $i++) {
  if($points[$i] - $start >= $min_block &&
     $ref_len - $points[$i] >= $min_block)
  {
    push(@points_tmp, $points[$i]);
    $start = $points[$i];
  }
}
@points = @points_tmp;

# Make blocks
$start = 0;
for(my $i = 0; $i < @points; $i++) {
  push(@blocks, make_block($start, $points[$i]));
  $start = $points[$i];
}

# Last block
if($start < $ref_len) {
  push(@blocks, make_block($start, $ref_len));
}

# Shuffle blocks
my @tmp_mix = shuffle(@blocks);
for my $b (@tmp_mix) { $b->{'strand'} = rand(100) < 50 ? REVERSE : FORWARD; }

# Merge adjacent blocks
my @mix = ($tmp_mix[0]);
for(my $i = 1; $i < @tmp_mix; $i++) {
  my ($b,$c) = ($tmp_mix[$i-1],$tmp_mix[$i]);
  my $same_strand = ($b->{'strand'} == $c->{'strand'});
  if($same_strand && $c->{'strand'} == FORWARD && $b->{'end'}+1 == $c->{'start'})
  {
    $mix[$#mix]->{'end'} = $c->{'end'};
  }
  elsif($same_strand && $c->{'strand'} == REVERSE && $b->{'start'} == $c->{'end'}+1)
  {
    $mix[$#mix]->{'start'} = $c->{'start'};
  }
  else {
    push(@mix, $c);
  }
}

# Save FASTA
print FA ">$title shuffled\n";
for my $b (@mix) {
  my $seq = substr($seq,  $b->{'start'}, $b->{'length'});
  if($b->{'strand'}) { $seq = rev_comp($seq); }
  print FA "$seq";
}
print FA "\n";

# Save block order
print TXT "# start  end  strand  (1-based coords)\n";
for my $m (@mix) {
  print TXT "".($m->{'start'}+1)."\t".
               ($m->{'end'}+1)."\t".
               ($m->{'strand'} == FORWARD ? '+' : '-')."\n";
}

close(FA);
close(TXT);

sub make_block
{
  my ($start,$next) = @_;
  my $end = $next - 1;
  my $length = $end - $start + 1;
  return {'start' => $start, 'end' => $end, 'length' => $length};
}
