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
"Usage: ./block_to_brkpnts.pl <ref.fa> <blocks.txt>\n" .
"  Convert blocks to breakpoints\n" .
"  Blocks should be (0-based): <start> tab <end> tab <strand=+->\n" .
"    e.g. 10\t100\t+\n";

  exit(-1);
}

if(@ARGV != 2) { print_usage(); }
my ($ref_path, $blocks_path) = @ARGV;

# Load sequence
my $fastn = open_fastn_file($ref_path);
my ($title,$ref) = $fastn->read_next();
if(!defined($title)) { die("Woops, can't read sequence: $ref_path"); }
close_fastn_file($fastn);

$ref = uc($ref);
($title) = ($title =~ /^(\S*)/);
my $ref_len = length($ref);
print STDERR "ref: '$title'\n";


open(BLOCKS, "$blocks_path") or die("Cannot open $blocks_path");
my $line;
my @mix = ();
while(defined($line = <BLOCKS>)) {
  chomp($line);
  if($line =~ /^#/) {}
  elsif($line =~ /^(\d+)\t(\d+)\t([\+\-])$/) {
    if($2 < $1) { die("Mixed up coords: $line"); }
    my ($start,$end,$s) = ($1,$2,$3);
    $start--; $end--; # convert to 0-based
    push(@mix, {'start' => $start, 'end' => $end,
                'strand' => ($s eq '+' ? FORWARD : REVERSE),
                'length' => $end - $start + 1});
  } else {
    die("Bad line: $line");
  }
}
close(BLOCKS);

my $max_shift = 0;
my $min_block = min(map {$_->{'length'}} @mix);
my $max_block = max(map {$_->{'length'}} @mix);

# Print 1-based coord file
for(my $i = 0; $i+1<@mix; $i++) {
  my $b = $mix[$i];
  my $c = $mix[$i+1];
  # Print in both directions
  print_breakpoint($b->{'strand'} ? $b->{'start'} : $b->{'end'},   $b->{'strand'},
                   $c->{'strand'} ? $c->{'end'}   : $c->{'start'}, $c->{'strand'});
  print "\t";
  print_breakpoint($c->{'strand'} ? $c->{'end'}   : $c->{'start'}, !$c->{'strand'},
                   $b->{'strand'} ? $b->{'start'} : $b->{'end'},   !$b->{'strand'});
  print "\n";
}

print STDERR "min block: $min_block\n";
print STDERR "max block: $max_block\n";
print STDERR "max shift: $max_shift\n";

if($max_shift > $min_block) {
  print STDERR "Fail: max shift > min block size\n";
  exit(-1);
}

# Shift breakpoints to the right due to matching sequence
# base $pos0:$strand0 is adjacent to $pos1:$strand1
sub print_breakpoint
{
  my ($pos0, $strand0, $pos1, $strand1) = @_;
  my $pos0_init =  $pos0;
  if($strand0 == FORWARD) {
    while($pos0+1 < $ref_len &&
          get_ref_base($pos0+1, $strand0) eq get_ref_base($pos1, $strand1)) {
      $pos0++;
      if($strand1 == FORWARD) { $pos1++; } else { $pos1--; }
    }
  } else {
    while($pos0 > 0 &&
          get_ref_base($pos0-1, $strand0) eq get_ref_base($pos1, $strand1)) {
      $pos0--;
      if($strand1 == FORWARD) { $pos1++; } else { $pos1--; }
    }
  }
  $max_shift = max($max_shift, abs($pos0-$pos0_init));
  # Print 1-based coords
  print "$title:".($pos0+1).":".($strand0?'-':'+')."\t".
        "$title:".($pos1+1).":".($strand1?'-':'+')."\t";
}

sub get_ref_base
{
  my ($pos, $strand) = @_;
  my $c = substr($ref, $pos, 1);
  return $strand == FORWARD ? $c : rev_comp($c);
}
