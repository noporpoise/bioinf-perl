#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first sum shuffle);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use GeneticsModule;
use UsefulModule;

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
close_fastn_file($fastn);
$seq = uc($seq);

if(!defined($title)) { die("Woops, can't read sequence: $ref_path"); }

my $ref_len = length($seq);

# Generate N breakpoints
my @points = sort {$a <=> $b} map {int(rand($ref_len-1))} 1..$N;

my @blocks = ();
my $start = 0;
my $min_block = 40;

for(my $i = 0; $i < @points && $start < $ref_len; $i++) {
  push(@blocks, make_block($start, $points[$i]));
  $start = $points[$i]+1;
}

# Last block
if($start < $ref_len) {
  push(@blocks, make_block($start, $ref_len-1));
}

# Shuffle blocks
# strand: 0 => FORWARD; 1 => REVERSE
my @mix = shuffle(@blocks);
for my $b (@mix) { $b->{'strand'} = rand(100) < 50 ? 1 : 0; }

# Save FASTA
print FA ">$title shuffled\n";
for my $b (@mix) {
  my $seq = substr($seq, $b->{'start'}, $b->{'len'});
  if($b->{'strand'}) { $seq = rev_comp($seq); }
  print FA "$seq";
}
print FA "\n";

# Save description
# for my $b (@blocks) {
#   print "$b->{'start'} $b->{'end'} $b->{'fw'}\n";
# }

# print breakpoint from block i to i+1
# for(my $i = 0; $i+1 < @mix; $i++) {
#   my $b = $mix[$i];
#   my $c = $mix[$i+1];
#   my $x = breakpoint_str($b->{'start'}+1, $b->{'end'}+1,    $b->{'strand'},
#                          $c->{'start'}+1, $c->{'end'}+1,    $c->{'strand'});
#   my $y = breakpoint_str($c->{'end'}+1,   $c->{'start'}+1, !$c->{'strand'},
#                          $b->{'end'}+1,   $b->{'start'}+1, !$b->{'strand'});
#   print TXT "$x or $y\n";
# }

# Print 1-based coord file
for my $b (@mix) {
  print TXT ($b->{'start'}+1)."\t".($b->{'end'}+1)."\t".($b->{'strand'}?'-':'+')."\n";
}

close(FA);
close(TXT);

sub breakpoint_str
{
  my ($start0,$end0,$rev0,$start1,$end1,$rev1) = @_;
  my $txt = ($rev0 ? ".*-$start0:-" : ".*-$end0:+") . " " .
            ($rev1 ? "$end1-.*:-" : "$start1-.*:-");
  return $txt;
}

sub make_block
{
  my ($start,$end) = @_;
  my $len = $end-$start+1;
  return {'start' => $start, 'end' => $end, 'len' => $len};
}
