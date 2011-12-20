#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(max);

use FASTNFile;
use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "usage: ./vcf_filter_by_ref.pl <K> <in.vcf> <ref files> ...\n";
  print STDERR "  Filter out variants that have Ns in the reference\n";
  print STDERR "  (useful for filtering out variants in repeats by using masked ref)\n";
  print STDERR "  If <in.vcf> is '-' then reads from STDIN\n";
  print STDERR "  <K> is flank size to check\n";
  exit;
}

if(@ARGV < 3)
{
  print_usage();
}

my $flank_size = shift;
my $vcf_file = shift;
my @ref_files = @ARGV;

if($flank_size !~ /^\d+$/) {
  print_usage("Flank size <k> must be a positive integer");
}

#
# Open VCF Handle
#
my $vcf_handle;

if($vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}


#
# Load reference files
#
my ($ref_genomes_hashref) = read_all_from_files(@ref_files);

my %ref_genomes = %$ref_genomes_hashref;

# Correct '1..' to 'chr1...' etc
#my ($key,$value)
while(my ($key,$value) = each(%ref_genomes))
{
  my $new_key = get_clean_chr_name($key);

  if($new_key ne $key)
  {
    $ref_genomes{$new_key} = uc($ref_genomes{$key});
    delete($ref_genomes{$key});
  }
  else {
    $ref_genomes{$key} = uc($ref_genomes{$key});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

print $vcf->get_header();

my $vcf_entry;
my %chrom_not_in_ref = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $chr = $vcf_entry->{'CHROM'};
  my $start = $vcf_entry->{'true_POS'} - $flank_size - 1;
  my $length = $flank_size + length($vcf_entry->{'true_REF'}) + $flank_size; 

  if($start < 0)
  {
    $length += $start; # reduce length
    $start = 0;
  }

  if(!defined($ref_genomes{$chr}))
  {
    if(!defined($chrom_not_in_ref{$chr}))
    {
      $chrom_not_in_ref{$chr} = 1;

      print STDERR "Warning: chromosome '$chr' not in reference - " .
                   "filtering out variants\n";
    }
  }
  elsif($start + $length <= length($ref_genomes{$chr}) &&
        substr($ref_genomes{$chr}, $start, $length) =~ /^[acgt]*$/i)
  {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
