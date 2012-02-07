#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_sum_svlens.pl [file.vcf]\n";
  print STDERR "  Sum polarised SVLENs\n";
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

my $num_of_variants = 0;
my $num_of_variants_polarised = 0;

my $num_mnps = 0;
my $sum_mnps = 0;

my $num_insertions = 0;
my $sum_insertions = 0;
my $num_deletions = 0;
my $sum_deletions = 0;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;
  my $svlen = $vcf_entry->{'INFO'}->{'AALEN'};
  
  if(defined($svlen))
  {
    $num_of_variants_polarised++;

    if($svlen > 0)
    {
      $num_insertions++;
      $sum_insertions += $svlen;
    }
    elsif($svlen < 0)
    {
      $num_deletions++;
      $sum_deletions += $svlen;
    }
    else
    {
      $num_mnps++;
      $sum_mnps += length($vcf_entry->{'true_REF'});
    }
  }
}

my $num_of_indels = $num_insertions + $num_deletions;

print "".num2str($num_of_variants_polarised)." / ".num2str($num_of_variants)." ".
      sprintf("%.2f", ($num_of_variants_polarised / $num_of_variants)) . "% " .
      " variants polarised\n";

print "Ins/Del per indel: ".(($sum_insertions+$sum_deletions) / $num_of_indels)."bp\n";

print "Total insertions: $sum_insertions for $num_insertions (".($sum_insertions/$num_insertions).")\n";
print "Total deletions: $sum_deletions for $num_deletions (".($sum_deletions/$num_deletions).")\n";
print "Total MNPs: $sum_mnps for $num_mnps (".($sum_mnps/$num_mnps).")\n";

close($vcf_handle);
