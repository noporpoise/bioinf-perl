#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_sort.pl [files..]
  Sort FASTA/Q file by read/chrom names.\n";

  exit;
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

  my %reads;
  my %read_quals;
  my ($title, $seq, $quals);

  while((($title, $seq, $quals) = $fastn->read_next()) && defined($title))
  {
    if(defined($reads{$title})) {
      print STDERR "fastn_sort.pl:Warning: duplicate reads with name '$title'\n";
    }
    $reads{$title} = $seq;
    $read_quals{$title} = $quals;
  }

  foreach my $title (sort keys %reads)
  {
    $fastn->print_entry($title, $reads{$title}, $read_quals{$title}, 80);
  }

  close_fastn_file($fastn);
}
