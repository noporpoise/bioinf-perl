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
  
  print STDERR "Usage: ./plot_debruijn.pl [--revcmp|--points] <kmer_size> <in.fa in.fq ...>\n";
  print STDERR "  If --revcmp is used, the lowest kmer will be used (2:1 mapping)\n";

  exit;
}

my $use_keys = 0;
my $use_points = 0;

while(@ARGV > 1) {
  if($ARGV[0] =~ /^-?-r(ev)?c?(mp)?$/i) { shift; $use_keys = 1; }
  elsif($ARGV[0] =~ /^-?-p(oints)?$/i) { shift; $use_points = 1; }
  else {last;}
}

if(@ARGV < 1) { print_usage(); }

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
print "  node [" . ($use_points ? "shape=point label=none" : "shape=plaintext") ."]\n";

if($use_keys) {
  print "  edge [dir=both arrowhead=none arrowtail=none]\n";
}

for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq);

  while((($title,$seq) = $fastn->read_next()) && defined($title))
  {
    for my $seq (split(/[^acgt]+/i, $seq)) {
      parse_read(uc($seq));
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
  my $prev_orient = 0;
  my $next_orient = 0;

  if($use_keys) { ($prev_kmer, $prev_orient) = get_kmer_key($prev_kmer); }

  for(my $i = 1; $i <= length($read)-$kmer_size; $i++)
  {
    my $next_kmer = substr($read, $i, $kmer_size);

    if($use_keys) {
      ($next_kmer, $next_orient) = get_kmer_key($next_kmer);
      print "$prev_kmer:".($prev_orient ? 'w' : 'e')." -> " .
            "$next_kmer:".($next_orient ? 'e' : 'w')."\n";
      ($prev_kmer,$prev_orient) = ($next_kmer,$next_orient);
    }
    else {
      print "$prev_kmer -> $next_kmer\n";
      $prev_kmer = $next_kmer;
    }
  }
}

sub get_kmer_key
{
  my $kmer = shift;
  my $rc = rev_comp($kmer);
  return ($kmer lt $rc ? ($kmer, 0) : ($rc, 1));
}
