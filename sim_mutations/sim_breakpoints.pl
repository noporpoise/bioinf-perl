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
my @points = sort {$a <=> $b} map {int(rand($ref_len-1))} 1..$N;

my @blocks = ();
my $start = 0;
my $min_block = 40;

# Remove blocks too short
print STDERR "Removing blocks shorter than $min_block\n";
my @points_tmp = ();
for(my $i = 0; $i < @points; $i++) {
  if($points[$i] - $start + 1 >= $min_block &&
     $ref_len - $points[$i] >= $min_block)
  {
    push(@points_tmp, $points[$i]);
    $start = $points[$i]+1;
  }
}
@points = @points_tmp;

# Make blocks
for(my $i = 0; $i < @points && $start < $ref_len; $i++) {
  push(@blocks, make_block($start, $points[$i]));
  $start = $points[$i]+1;
}

# Last block
if($start < $ref_len) {
  push(@blocks, make_block($start, $ref_len-1));
}

# Shuffle blocks
my @mix = shuffle(@blocks);
for my $b (@mix) { $b->{'strand'} = rand(100) < 50 ? REVERSE : FORWARD; }

# Save FASTA
print FA ">$title shuffled\n";
for my $b (@mix) {
  my $seq = substr($seq,  $b->{'start'}, $b->{'len'});
  if($b->{'strand'}) { $seq = rev_comp($seq); }
  print FA "$seq";
}
print FA "\n";

# Save description
# Print 1-based coord file
for(my $i = 0; $i+1<@mix; $i++) {
  my $b = $mix[$i];
  my $c = $mix[$i+1];
  # Print in both directions
  print_breakpoint($b->{'strand'} ? $b->{'start'} : $b->{'end'},   $b->{'strand'},
                   $c->{'strand'} ? $c->{'end'}   : $c->{'start'}, $c->{'strand'});
  print TXT "\t";
  print_breakpoint($c->{'strand'} ? $c->{'end'}   : $c->{'start'}, !$c->{'strand'},
                   $b->{'strand'} ? $b->{'start'} : $b->{'end'},   !$b->{'strand'});
  print TXT "\n";
}

close(FA);
close(TXT);

sub make_block
{
  my ($start,$end) = @_;
  my $len = $end-$start+1;
  return {'start' => $start, 'end' => $end, 'len' => $len};
}

# Shift breakpoints to the right due to matching sequence
sub print_breakpoint
{
  my ($pos0, $strand0, $pos1, $strand1) = @_;
  while($strand0 == FORWARD && $pos0+1 < $ref_len &&
        get_ref_base($pos0+1, $strand0) eq get_ref_base($pos1, $strand1)) {
    $pos0++;
    if($strand1 == FORWARD) { $pos1++; } else { $pos1--; }
  }
  while($strand0 == REVERSE && $pos0 > 0 &&
        get_ref_base($pos0-1, $strand0) eq get_ref_base($pos1, $strand1)) {
    $pos0--;
    if($strand1 == FORWARD) { $pos1++; } else { $pos1--; }
  }
  # Print 1-based coords
  print TXT "$title:".($pos0+1).":".($strand0?'-':'+')."\t".
            "$title:".($pos1+1).":".($strand1?'-':'+')."\t";
}

sub get_ref_base
{
  my ($pos, $strand) = @_;
  my $c = substr($seq, $pos, 1);
  return $strand == FORWARD ? $c : rev_comp($c);
}
