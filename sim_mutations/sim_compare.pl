#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

use GeneticsModule;
use UsefulModule;
use FASTNFile;
use VCFFile;

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

while(defined($entry = $vcf->read_entry()))
{
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if(!defined($vars{$vkey})) { $vars{$vkey} = {}; }
  map {$vars{$vkey}->{$_} = 1} @alleles;
}

close($fh);

#
# Parse result.vcf
#
my $num_true_positives = 0;
my $num_false_positives = 0;

open($fh, $result_vcf) or die("Cannot open result vcf: $result_vcf");
$vcf = new VCFFile($fh);

while(defined($entry = $vcf->read_entry()))
{
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if(!defined($vars{$vkey})) { $num_false_positives++; }
  else { $num_true_positives++; }
}

close($fh);

my $max_true_positives = scalar(keys(%vars));

print "Discovered: ".pretty_fraction($num_true_positives, $max_true_positives)."\n";
print "False positives: ".num2str($num_false_positives)."\n";

sub all { $_ || return 0 for @_; 1 }

sub resolve_var
{
  my ($lflank,$rflank,@alleles) = @_;
  my ($vkey,$rev) = get_var_key($lflank,$rflank);

  # Trim Padding left base
  my $base = substr($alleles[0],0,1);
  if(all(map {substr($_,0,1) eq $base} @alleles)) {
    @alleles = map {substr($_,1)} @alleles;
  }

  if($rev)
  {
    ($lflank,$rflank) = (rev_comp($rflank), rev_comp($lflank));
    @alleles = map {rev_comp($_)} @alleles;
  }

  return ($vkey,@alleles);
}

# Variant key is defined by <key1><key2>
# where <key1> is less than <key2>
# and the keys are the kmer keys of the start of the flanks
# if <key1> == <key2>, then the variant is `forward` when <kmer1>==<key1>
sub get_var_key
{
  my ($kmer1,$kmer2) = @_;
  my ($key1,$key2) = map {get_key($_)} @_;

  if($key2 gt $key1 || ($kmer1 eq $kmer2 && $kmer1 ne $key1))
  {
    # reverse variant
    return ($key2.$key1, 1);
  }
  return ($key1.$key2, 0);
}

sub get_key
{
  my $rev_comp = rev_comp($_[0]);
  return ($rev_comp lt $_[0] ? $rev_comp : $_[0]);
}
