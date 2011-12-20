#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw[min max];

sub printUsage
{
  if(@_) {
    print STDERR "Error: $_[0]\n";
  }
  
  print "usage: ./vcf_indel_sizes.pl [-p] <hist_size> [file]\n";
  print "  -p  only include indels that PASS all filters\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 3)
{
  printUsage();
}

my $filter = 0;

if($ARGV[0] eq "-p" || $ARGV[0] eq "-P")
{
  $filter = 1;
  shift;
}

my $hist_size = shift;
my $vcf_file = shift;

if(!defined($hist_size)) {
  printUsage();
}

my $handle;

if(defined($vcf_file)) {
  open($handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  printUsage("Must specify a file or pipe in VCF file");
}

my $line;

# Skip VCF headers
while(defined($line = <$handle>) && $line =~ /#/) {}

my @indel_dist = ();
my $num_of_vars = 0;

while(defined($line))
{
  my @cols = split(/\t/, $line);

  if(!$filter || $cols[6] =~ /PASS/)
  {
    my $ref_len = length($cols[3]);
  
    # Multiple alt alleles separated by commas
    my @alt_alleles = split(/,/, $cols[4]);
  
    for my $alt_allele (@alt_alleles)
    {
      my $alt_len = length($cols[4]);
      my $diff = $alt_len - $ref_len;

      if(min($ref_len, $alt_len) == 1)
      {
        $num_of_vars++;

        if(abs($diff) <= $hist_size)
        {
          $indel_dist[$diff + $hist_size]++;
        }
      }
    }
  }

  # Read next line
  $line = <$handle>;
}

print "num_of_vars,$num_of_vars\n";

for(my $i = -$hist_size; $i <= $hist_size; $i++)
{
  print "$i," .
       (defined($indel_dist[$i+$hist_size]) ? $indel_dist[$i+$hist_size] : '0') .
       "\n";
}

close($handle);

