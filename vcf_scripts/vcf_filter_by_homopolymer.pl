#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./vcf_filter_by_homopolymer.pl <max_run> [in.vcf]\n";
  print STDERR "  Filter out homopolymer runs LONGER than <max_run> (min = 1)\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $max_run = shift;
my $vcf_file = shift; # optional

if($max_run !~ /^\d+$/ || $max_run < 1) {
  print_usage("<max_run> must be >= 1");
}

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
print $vcf->get_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  if(length($vcf_entry->{'true_REF'}) != 0 &&
     length($vcf_entry->{'true_ALT'}) != 0)
  {
    die("Not a set of clean indels: use vcf_get_clean_indels.pl");
  }

  if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
     !defined($vcf_entry->{'INFO'}->{'right_flank'}))
  {
    die("'left_flank' or 'right_flank' INFO field missing: use vcf_add_flanks.pl");
  }

  my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
  my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

  my $insert = (length($vcf_entry->{'true_REF'}) > 0 ? $vcf_entry->{'true_REF'}
                                                     : $vcf_entry->{'true_ALT'});

  
  my $insert_first_char = substr($insert, 0, 1);
  my $all_one_base_insert = $insert_first_char x length($insert);

  my $all_bases_equal = ($insert =~ /$all_one_base_insert/i);

  my $matching_left_bases = 0;
  my $matching_right_bases = 0;

  if(reverse($left_flank) =~ /^($insert_first_char+)/) {
    $matching_left_bases = length($1);
  }
  
  if($right_flank =~ /^($insert_first_char+)/) {
    $matching_right_bases = length($1);
  }

  my $hrun = $matching_left_bases + $matching_right_bases + length($insert);

  #print "bases_equal=$all_bases_equal, " .
  #      "matching_left_bases=$matching_left_bases," .
  #      "matching_right_bases=$matching_right_bases, hrun=$hrun\n";

  if(!$all_bases_equal || $hrun <= $max_run) {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
