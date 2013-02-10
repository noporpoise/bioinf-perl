#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use GeneticsModule;

#
# Read in fasta/fastq files, produce de Bruijn graph in dot format
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./plot_debruijn.pl [--revcmp] <kmer_size> <in.fa in.fq ...>\n";
  print STDERR "  If --revcmp is used, the lowest kmer will be used (2:1 mapping)\n";

  exit;
}

my $use_keys = 0;
if(@ARGV < 1) { print_usage(); }
if($ARGV[0] =~ /^-?-r(ev)?c?(mp)?$/i) { shift; $use_keys = 1; }

my $kmer_size = shift;
my @files = @ARGV;

if($kmer_size !~ /^\d+$/ || $kmer_size < 1)
{
  print_usage("Invalid kmer_size value '$kmer_size'");
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

  if(length($read) < $kmer_size)
  {
    return;
  }

  my $prev_kmer = substr($read, 0, $kmer_size);
  if($use_keys) { $prev_kmer = rev_comp_key($prev_kmer); }

  for(my $i = 1; $i <= length($read)-$kmer_size; $i++)
  {
    my $new_kmer = substr($read, $i, $kmer_size);
    if($use_keys) { $new_kmer = rev_comp_key($new_kmer); }
    print "$prev_kmer -> $new_kmer\n";
    $prev_kmer = $new_kmer;
  }
}

sub rev_comp_key
{
  my $kmer = shift;
  my $rc = rev_comp($kmer);
  return $kmer lt $rc ? $kmer : $rc;
}
