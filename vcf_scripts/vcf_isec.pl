#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule; # num2str
use VCFFile;
use RefGenome;

sub print_usage
{
  for my $err (@_) { print STDERR "Error: $err\n"; }
  print STDERR "Usage: ./vcf_isec.pl <in.vcf> <in2.vcf> ..\n";
  exit;
}

if(@ARGV < 2) { print_usage(); }
my @files = @ARGV;

#
# Open VCF Handles
#
my @handles = ();
my @vcfs = ();

for my $file (@files) {
  my $vcf_handle;
  if(defined($file) && $file ne "-")
  {
    open($vcf_handle, $file)
      or print_usage("Cannot open VCF file '$file'\n");
  }
  elsif(-p STDIN) {
    # STDIN is connected to a pipe
    open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
  }
  else {
    print_usage("Must specify or pipe in a VCF file");
  }
  push(@handles, $vcf_handle);
  push(@vcfs, new VCFFile($vcf_handle));
}

# First vcf file is the one we take entries from,
# then check that they are in the others files
my $first_vcf = shift(@vcfs);
$first_vcf->print_header();

my $nfiles = @vcfs;
# my @entries = ([] x $nfiles);
my ($nprinted, $nentries) = (0, 0);

my @entries;
for(my $i = 0; $i < $nfiles; $i++) { $entries[$i] = []; }

while(defined(my $vcf_entry = $first_vcf->read_entry()))
{
  # $first_vcf->print_entry($vcf_entry);
  my $i; $nentries++;
  for($i = 0; $i < $nfiles; $i++) {
    read_vcf_entries($vcfs[$i], $entries[$i], $vcf_entry->{'CHROM'}, $vcf_entry->{'POS'});
  }
  for($i = 0; $i < $nfiles && entry_in_list($vcf_entry, $entries[$i]); $i++) {}
  if($i == $nfiles) { $nprinted++; $first_vcf->print_entry($vcf_entry); }
}

print "".pretty_fraction($nprinted, $nentries)." entries printed\n";

close($first_vcf->{'_handle'});
for my $vcf (@vcfs) { close($vcf->{'_handle'}); }

sub entry_in_list
{
  my ($entry,$arr) = @_;
  for(my $i = 0; $i < @$arr; $i++) {
    if($arr->[$i]->{'CHROM'} eq $entry->{'CHROM'} &&
       $arr->[$i]->{'POS'} eq $entry->{'POS'} &&
       $arr->[$i]->{'true_REF'} eq $entry->{'true_REF'}) {
      return 1;
    }
  }
  return 0;
}

sub read_vcf_entries
{
  my ($vcf,$arr,$chr,$pos) = @_;

  # Remove old entries
  while(scalar(@$arr) > 0 &&
       ($arr->[0]->{'CHROM'} lt $chr ||
        ($arr->[0]->{'CHROM'} eq $chr && $arr->[0]->{'POS'} < $pos)))
  {
    shift(@$arr);
  }

  # Get new ones
  my $last = @$arr == 0 ? undef : $arr->[@$arr-1];
  if(!defined($last) || ($last->{'CHROM'} eq $chr && $last->{'POS'} == $pos))
  {
    while(defined(my $vcf_entry = $vcf->read_entry()))
    {
      if($vcf_entry->{'CHROM'} ge $chr && $vcf_entry->{'POS'} >= $pos) {
        push(@$arr, $vcf_entry);
        if($vcf_entry->{'CHROM'} gt $chr || $vcf_entry->{'POS'} > $pos) { last; }
      }
    }
  }
}
