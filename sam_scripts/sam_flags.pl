#!/usr/bin/perl

use strict;
use warnings;

my @flagtxt =
 ("  0x1: Template has multiple fragments (mate paired)",
  "  0x2: Each frament properly aligned",
  "  0x4: Fragment unmapped",
  "  0x8: Next fragment in the template unmapped",
  " 0x10: SEQ being reverse complemented",
  " 0x20: SEQ of the next fragment in the template being reversed",
  " 0x40: The first fragment in the template (first in pair)",
  " 0x80: The last fragment in the template (second in pair)",
  "0x100: Secondary alignment",
  "0x200: Not passing quality controls",
  "0x400: PCR or optical duplicate");

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "usage: ./sam_flags.pl <number>\n";
  print STDERR join('', map {" $_\n"} @flagtxt);
  exit(-1);
}

if(@ARGV != 1) { print_usage();}

my $flags = shift;

if($flags & 0x1)   { print $flagtxt[0]."\n"; }
if($flags & 0x2)   { print $flagtxt[1]."\n"; }
if($flags & 0x4)   { print $flagtxt[2]."\n"; }
if($flags & 0x8)   { print $flagtxt[3]."\n"; }
if($flags & 0x10)  { print $flagtxt[4]."\n"; }
if($flags & 0x20)  { print $flagtxt[5]."\n"; }
if($flags & 0x40)  { print $flagtxt[6]."\n"; }
if($flags & 0x80)  { print $flagtxt[7]."\n"; }
if($flags & 0x100) { print $flagtxt[8]."\n"; }
if($flags & 0x200) { print $flagtxt[9]."\n"; }
if($flags & 0x400) { print $flagtxt[10]."\n";}
