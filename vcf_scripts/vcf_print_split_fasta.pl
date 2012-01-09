#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_split_fasta.pl [file.vcf]\n";
  print STDERR "  Print FASTA of 5' branch, ancestral, alt allele & 3' branch\n";
  print STDERR "  On different lines.  Needs AA tag to be 0 or 1\n";
  exit;
}

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
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

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Works for multiple ALT frequencies separated by commas
  my $id = $vcf_entry->{'ID'};
  my $ref = $vcf_entry->{'true_REF'};
  my $alt = $vcf_entry->{'true_ALT'};

  my $left = $vcf_entry->{'INFO'}->{'left_flank'};
  my $right = $vcf_entry->{'INFO'}->{'left_flank'};

  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  if($aa eq "0" || $aa eq "1")
  {
    print ">".$id."_5p\n";
    print "$left\n";
    print ">".$id."_anc\n";
    print "".($aa eq "0" ? $ref : $alt) . "\n";
    print ">".$id."_alt\n";
    print "".($aa eq "1" ? $ref : $alt) . "\n";
    print ">".$id."_3p\n";
    print "$right\n";
  }
}

close($vcf_handle);
