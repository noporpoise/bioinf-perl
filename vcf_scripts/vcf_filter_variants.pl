#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use UsefulModule; # for num2str

sub print_usage
{
  for my $err (@_)
  {
    chomp($err);
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_filter_variants.pl [--invert] <INDEL|CLEAN_INDEL|SNP|MNP|NP> [in.vcf]\n";
  print STDERR "  Prints variants were SVLEN != 0\n";
  print STDERR "  --invert      Everything but option\n";
  print STDERR "  \n";
  print STDERR "  Options:\n";
  print STDERR "   INDEL        SVLEN != 0\n";
  print STDERR "   CLEAN_INDEL  (SVLEN != 0) & (one allele 0bp)\n";
  print STDERR "   SNP          both alleles 1bp\n";
  print STDERR "   MNP          (SVLEN == 0) & (allele length > 1)\n";
  print STDERR "   NP           SVLEN == 0\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 3)
{
  print_usage();
}

my $filter_invert = 0;

if(@ARGV > 1 && $ARGV[0] =~ /^-?-c(lean)?$/i)
{
  $filter_invert = 1;
  shift;
}

my $filter = shift;
$filter = uc($filter);

my @filters = qw(INDEL CLEAN_INDEL SNP MNP NP);

if(!grep {$_ eq $filter} @filters)
{
  print_usage("Unknown option: '$filter' " .
              "(expected one of: ".join(",",@filters).")");
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
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
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $print;

  if($filter eq "INDEL")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0);
  }
  elsif($filter eq "CLEAN_INDEL")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0 &&
              (length($vcf_entry->{'true_REF'}) == 0 ||
               length($vcf_entry->{'true_ALT'}) == 0));
  }
  elsif($filter eq "SNP")
  {
    $print = (length($vcf_entry->{'true_REF'}) == 1 &&
              length($vcf_entry->{'true_ALT'}) == 1);
  }
  elsif($filter eq "MNP")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0 &&
              length($vcf_entry->{'true_REF'}) > 1);
  }
  elsif($filter eq "NP")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0);
  }

  if($print != $filter_invert)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_printed++;
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_of_printed / $num_of_variants;

print STDERR "vcf_filter_variants.pl: (".($filter_invert ? "not " : "")."$filter) " .
             num2str($num_of_printed) . " / " . num2str($num_of_variants) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
