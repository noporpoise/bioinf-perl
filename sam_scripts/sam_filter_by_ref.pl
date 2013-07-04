#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use RefGenome;
use UsefulModule;

sub print_usage
{
  foreach my $err (@_) {print STDERR "Error: $err\n";}

  print STDERR "" .
"Usage: ./sam_filter_by_ref.pl <in.sam> <ref1.fa ..>
  Filter out entries with N in the ref\n";
  exit(-1);
}

if(@ARGV < 1) { print_usage(); }

my $sam_file = shift;
my @ref_files = @ARGV;

if(@ref_files == 0) { print_usage("No references files given"); }

#
# Open VCF Handle
#
my $sam_handle;

if($sam_file ne "-") {
  open($sam_handle, $sam_file) or die("Cannot open VCF file '$sam_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($sam_handle, "<&=STDIN") or die("Cannot read pipe");
}
else {
  print_usage("Must specify or pipe in a SAM file");
}

#
# Load reference files
#
my $genome = new RefGenome(uppercase => 1);
$genome->load_from_files(@ref_files);

#
# Read SAM
#
my $line;
my ($total_num_entries, $num_of_printed_entries) = (0,0);
my %chrom_not_in_ref = ();
my %coords_out_of_bounds = ();

while(defined($line = <$sam_handle>))
{
  if($line =~ /^@/ || $line =~ /^\s*$/) { print $line; next; }
  chomp($line);

  my ($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext,
      $tlen, $seq, $qual, $extra) = split("\t", $line);

  if(!defined($qual)) { die("Invalid entry: $line"); }

  $total_num_entries++;

  if(!$genome->chr_exists($rname))
  {
    if(!defined($chrom_not_in_ref{$rname})) {
      $chrom_not_in_ref{$rname} = 1;
      print STDERR "sam_filter_by_ref.pl - Warning: chromosome '$rname' not " .
                   "in reference, filtering out entries\n";
    }
    next;
  }

  my $start = $pos-1;
  my $length = length($seq);
  my $chrom_length = $genome->get_chr_length($rname);

  if($start + $length > $chrom_length)
  {
    if(!defined($coords_out_of_bounds{$rname})) {
      $coords_out_of_bounds{$rname} = 1;
      print STDERR "sam_filter_by_ref.pl - Warning: entry '$qname' " .
                   "[$rname:$pos:$length] " .
                   "out of bounds of ref " .
                   "'".$genome->guess_chrom_fasta_name($rname)."' " .
                   "[length:$chrom_length]\n";
    }
    next;
  }

  my $ref_seq = $genome->get_chr_substr0($rname, $start, $length);

  if($ref_seq =~ /^[ACGT]*$/)
  {
    $num_of_printed_entries++;
    print $line;
  }
}

print STDERR "sam_filter_by_ref.pl: " .
             pretty_fraction($num_of_printed_entries, $total_num_entries) .
             " entries printed\n";

close($sam_handle);
