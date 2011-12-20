#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_allele_freq.pl [file.vcf]\n";
  print STDERR "  Add (or correct) allele frequency INFO tag in a VCF file\n";
  print STDERR "  Uses AF (Alt allele Frequency) tag: 0..1\n";
  exit;
}

if(@ARGV > 1) {
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

my @sample_names = $vcf->get_list_of_sample_names();

if(@sample_names == 0)
{
  print_usage("VCF doesn't have any genotyped samples in it");
}

my $additional_header = '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele ' .
                        'Frequency, for each ALT allele, in the same order as listed">';

print vcf_add_to_header($vcf->get_header(), $additional_header);

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Works for multiple ALT frequencies separated by commas
  my @alleles = split(/,/, $vcf_entry->{'ALT'});

  my %allele_counts = ();

  for my $sample (@sample_names)
  {
    if($vcf_entry->{$sample} =~ /^(\d+)\/(\d+)/)
    {
      $allele_counts{$1}++;
      $allele_counts{$2}++;
    }
    else {
      die("Invalid entry '$vcf_entry->{$sample}'");
    }
  }

  my @alt_freqs = ();
  my $total_ploidy = 2 * @sample_names;

  for(my $i = 1; $i <= @alleles; $i++)
  {
    my $allele_count = defined($allele_counts{$i}) ? $allele_counts{$i} : 0;
    push(@alt_freqs, $allele_count / $total_ploidy);
  }

  $vcf_entry->{'INFO'}->{'AF'} = join(",", @alt_freqs);
  
  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);
