#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;

#
# Read in fasta/fastq files, produce de Bruijn graph in dot format
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./plot_debruijn.pl <kmer> <in.fa in.fq ...>\n";

  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $kmer = shift;
my @files = @ARGV;

if($kmer !~ /^\d+$/ || $kmer < 1)
{
  print_usage("Invalid kmer value '$kmer'");
}

if(@files == 0)
{
  push(@files, "-");
}

print "digraph G {\n";

for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq);

  while((($title,$seq) = $fastn->read_next()) && defined($title))
  {
    for my $seq (split(/[^acgt]+/i, $seq)) {
      parse_read($seq);
    }
  }

  close_fastn_file($fastn);
}

print "}\n";

sub parse_read
{
  my ($read) = @_;

  if(length($read) < $kmer)
  {
    return;
  }

  my $prev_kmer = substr($read, 0, $kmer);

  for(my $i = 1; $i <= length($read)-$kmer; $i++)
  {
    my $new_kmer = substr($read, $i, $kmer);
    print "$prev_kmer -> $new_kmer\n";
    $prev_kmer = $new_kmer;
  }
}
