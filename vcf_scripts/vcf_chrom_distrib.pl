#!/usr/bin/perl

use strict;
use warnings;

use POSIX qw(ceil);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_chrom_distrib.pl <chrom_sizes.csv> <bin_size> [vcf]\n";
  print STDERR "  Get distribution of variants across the genome\n";
  print STDERR "  <chrom_sizes.csv> is a csv file of chromosome lengths\n";
  print STDERR "  <bin_size> is the number of kbp in each bin\n";
  print STDERR "  (reads from STDIN if [vcf] is omitted or is '-').\n";
  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 2 || @ARGV > 3)
{
  print_usage();
}

my $chrom_sizes_file = shift;
my $bin_size_kbp = shift;
my $vcf_file = shift;

if($bin_size_kbp !~ /^\d*\.?\d+$/)
{
  print_usage("Invalid bin size argument '$bin_size_kbp' - should be +ve kbp");
}

my $bin_size = $bin_size_kbp * 1000;

#
# Load chromosome lengths
#
my %chrom_lengths = ();

open(CHROMS, $chrom_sizes_file)
  or die("Cannot open chromosome length file '$chrom_sizes_file'");

my $chrom_line = <CHROMS>;

if($chrom_line !~ /^\w+,\d+$/)
{
  # May have read header line
  $chrom_line = <CHROMS>
}

while(defined($chrom_line))
{
  my ($chr, $length) = ($chrom_line =~ /^(\w+)[,\s]+(\d+)$/);
  $chr = get_clean_chr_name($chr);
  $chrom_lengths{$chr} = $length;

  $chrom_line = <CHROMS>;
}

close(CHROMS);

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my %bins_by_chr = ();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $chr = $vcf_entry->{'CHROM'};

  if(!defined($bins_by_chr{$chr}))
  {
    my $num_of_bins = ceil($chrom_lengths{$chr} / $bin_size);
    # Create an array of length $num_of_bins, full of zeros
    $bins_by_chr{$chr} = \(map {0} 1..$num_of_bins);
  }

  my $bin_index = int($vcf_entry->{'true_POS'} / $bin_size);

  $bins_by_chr{$chr}->[$bin_index]++;
}

$vcf->vcf_close();

print "chr,bin_index,bin_start,count\n";

for my $chr (sort {$a cmp $b} keys %bins_by_chr)
{
  my $bin_arr = $bins_by_chr{$chr};

  for(my $bin_index = 0; $bin_index < @$bin_arr; $bin_index++)
  {
    print "$chr,$bin_index," . ($bin_index * $bin_size) . "," .
          $bin_arr->[$bin_index] . "\n";
  }
}
