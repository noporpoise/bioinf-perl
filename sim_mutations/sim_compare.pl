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

  if(!defined($vars{$vkey})) {
    $vars{$vkey} = {'seen'=>0, 'alleles'=>\@alleles};
  }
  else {
  #   my $arr = $vars{$vkey}->{'alleles'};
    push(@{$vars{$vkey}->{'alleles'}}, @alleles);
  }
}
while(defined($entry = $vcf->read_entry()));

close($fh);

#
# Parse result.vcf
#
my $num_true_positives = 0;
my $num_false_positives = 0;
my $num_dupes = 0;

my $vars_matching_alleles = 0;
my $vars_missing_alleles = 0;
my $vars_false_alleles = 0;

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

    my ($num_matching, $num_missing, $num_false)
      = compare_alleles($vars{$vkey}->{'alleles'}, \@alleles);

    if($num_matching > 0) { $vars_matching_alleles++; }
    if($num_missing > 0) { $vars_missing_alleles++; }
    if($num_false > 0) { $vars_false_alleles++; }
  }
}

close($fh);

my $max_true_positives = scalar(keys(%vars));

print "Discovered: " . pretty_fraction($num_true_positives,
                                       $max_true_positives) . "\n";
print "False positives: " . num2str($num_false_positives) . "\n";
print "Dupes: " . num2str($num_dupes) . "\n";

print "Vars with at least one matching alleles: " .
      pretty_fraction($vars_matching_alleles, $num_true_positives) . "\n";
print "Vars missing alleles: " .
      pretty_fraction($vars_missing_alleles, $num_true_positives) . "\n";
print "Vars with extra alleles: " .
      pretty_fraction($vars_false_alleles, $num_true_positives) . "\n";

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
  ($lflank,$rflank,@alleles) = normalise_variant($kmer_size, $lflank, $rflank,
                                                 @alleles);

  my ($key0, $key1) = map {get_key($_)} (substr($lflank,-$kmer_size),
                                         substr($rflank,0,$kmer_size));

  return ($key0.$key1, @alleles);
}

# returns: (num alleles matching, num alleles missing, num false alleles)
sub compare_alleles
{
  my ($trutharr,$testarr) = @_;

  my %truthhsh = ();
  map {$truthhsh{$_} = 1} @$trutharr;

  my ($num_matched, $num_false) = (0,0);

  for my $a (@$testarr) {
    if(defined($truthhsh{$a})) { $num_matched++; }
    else { $num_false++; }
  }

  my $num_missing = scalar(@$trutharr) - $num_matched;

  # print "truth: ".join(',',@$trutharr)."\n";
  # print "test : ".join(',',@$testarr)."\n";
  # print "  (".join(',',($num_matched, $num_missing, $num_false)).")\n";

  return ($num_matched, $num_missing, $num_false);
}

