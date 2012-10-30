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
  
  print STDERR "" .
"Usage: ./resolve_debruijn.pl <kmer> <in.fa in.fq ...>
  
  Example:
    ./resolve_debruijn.pl 5 graphs/test.fa > graphs/test.dot
    dot -Tpng graphs/test.dot > graphs/test.png\n";

  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $kmer_size = shift;
my @files = @ARGV;

if($kmer_size !~ /^\d+$/ || $kmer_size < 1)
{
  print_usage("Invalid kmer value '$kmer_size'");
}

if(@files == 0)
{
  push(@files, "-");
}

# $graph{kmers}->{nextbases} = 1
my %graph;

# Build graph
for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq,$qual);

  while((($title,$seq,$qual) = $fastn->read_next()) && defined($title))
  {
    load_read($seq);
  }

  close_fastn_file($fastn);
}

# Resolve structure
for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq,$qual);

  while((($title,$seq,$qual) = $fastn->read_next()) && defined($title))
  {
    resolve_read($seq);
  }

  close_fastn_file($fastn);
}

#dump_graph();

sub load_read
{
  my ($read) = @_;

  if(length($read) < $kmer_size)
  {
    return;
  }

  for(my $i = 0; $i < length($read)-$kmer_size; $i++)
  {
    my $new_kmer = substr($read, $i, $kmer_size);
    my $next_base = substr($read, $i+$kmer_size, 1);

    if(!defined($graph{$new_kmer})) {
      $graph{$new_kmer} = {};
    }

    $graph{$new_kmer}->{$next_base}++;
    #print "$new_kmer [$next_base]\n";
  }
}

sub dump_graph
{
  print "digraph G {\n";
  for my $tmp_kmer (sort keys %graph)
  {
    my $next = substr($tmp_kmer, 1);
    for my $base (sort keys %{$graph{$tmp_kmer}})
    {
      print "$tmp_kmer -> $next$base [penwidth=".$graph{$tmp_kmer}->{$base}."];\n";
    }
  }
  print "}\n";
}

sub resolve_read
{
  my ($read) = @_;
  # 1) split by recurring kmer
  # 2) 
  print "read: $read\n";
  my ($khash, $karr) = get_hash_of_kmers($read);

  # Loop over kmers
  for(my $i = 0; $i < @$karr; $i++)
  {
    if($khash->{$karr->[$i]} > 1)
    {
      print "$i $karr->[$i] $khash->{$karr->[$i]}\n";
    }
  }
}

sub get_hash_of_kmers
{
  my ($read) = @_;
  my %kmers_hash = ();
  my @kmers = map {substr($read, $_, $kmer_size)} 0..(length($read)-$kmer_size-1);
  map {$kmers_hash{$_}++} @kmers;
  return (\%kmers_hash, \@kmers);
}
