#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule;
use VCFFile;
use SimLib;

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./sim_compare.pl <truth.vcf> <result.vcf>\n";
  exit(-1);
}

if(@ARGV != 2) { print_usage(); }
my $truth_vcf = shift;
my $result_vcf = shift;

# key => array of alleles
my %variants = ();

#
# Load truth VCF
#
my $fh;
open($fh, $truth_vcf) or die("Cannot open truth vcf: $truth_vcf");
my $vcf = new VCFFile($fh);
my $entry;
my %vars = ();

if(!defined($entry = $vcf->read_entry())) { die("Empty truth VCF"); }
my $kmer_size = length($entry->{'INFO'}->{'LF'});

do
{
  # print "a: ".$entry->{'ID'}."\n";
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if(!defined($vars{$vkey})) { $vars{$vkey} = {'seen'=>0}; }
}
while(defined($entry = $vcf->read_entry()));

close($fh);

#
# Parse result.vcf
#
my $num_true_positives = 0;
my $num_false_positives = 0;
my $num_dupes = 0;

open($fh, $result_vcf) or die("Cannot open result vcf: $result_vcf");
$vcf = new VCFFile($fh);

while(defined($entry = $vcf->read_entry()))
{
  # print "b: ".$entry->{'ID'}."\n";
  # if($entry->{'ID'} eq "var931") {last;}
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if(!defined($vars{$vkey})) { $num_false_positives++; }
  else {
    if($vars{$vkey}->{'seen'} == 0) { $num_true_positives++; }
    else { $num_dupes++; }
    $vars{$vkey}->{'seen'}++;
  }
}

close($fh);

my $max_true_positives = scalar(keys(%vars));

print "Discovered: " . pretty_fraction($num_true_positives,
                                       $max_true_positives) . "\n";
print "False positives: " . num2str($num_false_positives) . "\n";
print "Dupes: " . num2str($num_dupes) . "\n";

sub all { $_ || return 0 for @_; 1 }

sub resolve_var
{
  my ($lflank,$rflank,@alleles) = @_;

  # Remove padding base
  my $base = substr($alleles[0],0,1);
  if(all(map {substr($_,0,1) eq $base} @alleles)) {
    @alleles = map {substr($_,1)} @alleles;
  }

  # print_variant($lflank, $rflank, @alleles);
  ($lflank,$rflank,@alleles) = normalise_variant($kmer_size,$lflank,$rflank,@alleles);

  my ($key0, $key1) = map {get_key($_)} (substr($lflank,-$kmer_size),
                                         substr($rflank,0,$kmer_size));

  return ($key0.$key1, @alleles);
}



