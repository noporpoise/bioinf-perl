#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use FASTNFile;
use GeneticsModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_revcmp.pl [-w <wrap>] [files..]
  Reverse complement sequences.\n";

  exit;
}

my $linewrap = 80;
if(@ARGV > 2 && $ARGV[0] eq "-w") {
  shift;
  $linewrap = shift;
  if($linewrap !~ /^\d+$/) { print_usage("Invalid linewrap argument (-w)"); }
}

if(@ARGV == 0) { print_usage();}

my @files = @ARGV;

if(scalar(grep {$_ eq "-"} @files) > 1)
{
  print STDERR "Warning: reading from stdin more than once (multiple '-'s)\n";
}

for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title, $seq, $quals);

  while((($title, $seq, undef, $quals) = $fastn->read_next()) && defined($title))
  {
    $fastn->print_entry($title, rev_comp($seq), reverse($quals), $linewrap);
  }

  close_fastn_file($fastn);
}
