#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(first);
use File::Path qw(make_path);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule;
use GeneticsModule;
use VCFFile;
use FASTNFile;
use SimLib;
use IndexedString;

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./sim_compare.pl <truth.vcf> <result.vcf> " .
                 "<truth.out.vcf> <label> <false_positives.vcf> [genomes.fa]\n";
  exit(-1);
}

if(@ARGV < 5) { print_usage(); }
my $truth_vcf = shift;
my $result_vcf = shift;
my $out_vcf = shift;
my $label = shift;
my $failures_vcf = shift;
my @genome_files = @ARGV;

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
my $kmer_size;

my ($true_snp, $true_ins, $true_del, $true_inv) = (0,0,0,0);

if(defined($entry = $vcf->read_entry()))
{
  $kmer_size = length($entry->{'INFO'}->{'LF'});

  do
  {
    # print "a: ".$entry->{'ID'}."\n";
    my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                       $entry->{'INFO'}->{'RF'},
                                       split(',', $entry->{'ALT'}));

    if(!defined($vars{$vkey})) {
      $vars{$vkey} = {'seen'=>0, 'alleles'=>[]};
    }

    push(@{$vars{$vkey}->{'alleles'}}, @alleles);

    $vars{$vkey}->{'SNP'} += $entry->{'INFO'}->{'SNP'};
    $vars{$vkey}->{'INS'} += $entry->{'INFO'}->{'INS'};
    $vars{$vkey}->{'DEL'} += $entry->{'INFO'}->{'DEL'};
    $vars{$vkey}->{'INV'} += $entry->{'INFO'}->{'INV'};

    $true_snp += $vars{$vkey}->{'SNP'};
    $true_ins += $vars{$vkey}->{'INS'};
    $true_del += $vars{$vkey}->{'DEL'};
    $true_inv += $vars{$vkey}->{'INV'};
  }
  while(defined($entry = $vcf->read_entry()));
}
else { warn("Emtpy truth VCF"); }

close($fh);

#
# Parse result.vcf
#
my ($num_true_positives, $num_false_positives, $num_dupes) = (0,0,0);
my ($vars_matching_alleles, $vars_missing_alleles, $vars_false_alleles) = (0,0,0);
my ($num_snp_found, $num_ins_found, $num_del_found, $num_inv_found) = (0,0,0,0);

open($fh, $result_vcf) or die("Cannot open result vcf: $result_vcf");
$vcf = new VCFFile($fh);

while(defined($entry = $vcf->read_entry()))
{
  my $kmer_size2 = length($entry->{'INFO'}->{'LF'});
  if(defined($kmer_size) && $kmer_size != $kmer_size2) {
    # die("kmer sizes don't match [$kmer_size != $kmer_size2]");
  }
  $kmer_size = $kmer_size2;

  # print "b: ".$entry->{'ID'}."\n";
  # if($entry->{'ID'} eq "var931") {last;}
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if(defined($vars{$vkey}))
  {
    if($vars{$vkey}->{'seen'} == 0) { $num_true_positives++; }
    else { $num_dupes++; }

    my ($num_matching, $num_missing, $num_false)
      = compare_alleles($vars{$vkey}->{'alleles'}, \@alleles);

    if($num_matching > 0) { $vars_matching_alleles++; }
    if($num_missing > 0) {
      $vars_missing_alleles++;
      $vars{$vkey}->{'missing_alleles'} = 1;
    }
    if($num_false > 0) {
      $vars_false_alleles++;
      $vars{$vkey}->{'false_alleles'} = 1;
    }

    if($num_matching > 0) {
      $num_snp_found += $vars{$vkey}->{'SNP'};
      $num_ins_found += $vars{$vkey}->{'INS'};
      $num_del_found += $vars{$vkey}->{'DEL'};
      $num_inv_found += $vars{$vkey}->{'INV'};
    }

    $vars{$vkey}->{'seen'}++;
  }
  else { $num_false_positives++; }
}

close($fh);

my $max_true_positives = scalar(keys(%vars));

# 'False positive' is perhaps not the correct term for unexpected bubbles
print "Discovered: " . pretty_fraction($num_true_positives,
                                       $max_true_positives) . "\n";
print "Unexpected bubbles: " . num2str($num_false_positives) . "\n";
print "Dupes: " . num2str($num_dupes) . "\n";

print "---- breakdown ----\n";
print "Vars with at least one matching alleles: " .
      pretty_fraction($vars_matching_alleles, $num_true_positives) . "\n";
print "Vars missing alleles: " .
      pretty_fraction($vars_missing_alleles, $num_true_positives) . "\n";
print "Vars with extra alleles: " .
      pretty_fraction($vars_false_alleles, $num_true_positives) . "\n";

print "SNPs: " . pretty_fraction($num_snp_found, $true_snp) . "\n";
print "INSs: " . pretty_fraction($num_ins_found, $true_ins) . "\n";
print "DELs: " . pretty_fraction($num_del_found, $true_del) . "\n";
print "INVs: " . pretty_fraction($num_inv_found, $true_inv) . "\n";

#
# Print labelled truth vcf
#
my $out;
open($out, ">$out_vcf") or die("Cannot open failures VCF: $out_vcf");
open($fh, $truth_vcf) or die("Cannot open truth vcf: $truth_vcf");
$vcf = new VCFFile($fh);

while(defined($entry = $vcf->read_entry()))
{
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  if($vars{$vkey}->{'seen'} > 0) {
    $entry->{'INFO_flags'}->{$label} = 1;
  }
  $vcf->print_entry($entry, $out);
}

close($fh);
close($out);

# Load genomes to check alleles
my $search_genomes = new IndexedString($kmer_size);

for my $genome_file (@genome_files)
{
  my $fastn = open_fastn_file($genome_file);
  my ($title,$seq);
  while((($title,$seq) = $fastn->read_next()) && defined($title)) {
    $seq =~ s/[^ACGT]//gi;
    $search_genomes->add_strings($seq);
  }
  close_fastn_file($fastn);
}

#
# Print false positives vcf
#
open($out, ">$failures_vcf") or die("Cannot open failures VCF: $failures_vcf");
open($fh, $result_vcf) or die("Cannot open results output VCF: $result_vcf");
$vcf = new VCFFile($fh);

my $num_of_contigs = 0;
my $num_of_false_contigs = 0;

while(defined($entry = $vcf->read_entry()))
{
  my ($vkey, @alleles) = resolve_var($entry->{'INFO'}->{'LF'},
                                     $entry->{'INFO'}->{'RF'},
                                     split(',', $entry->{'ALT'}));

  my $lf = $entry->{'INFO'}->{'LF'};
  my $rf = $entry->{'INFO'}->{'RF'};

  my @contig_status = map {contig_exists($lf.$_.$rf) ? 1 : 0} @alleles;
  $entry->{'INFO'}->{'CONTIGS'} = join(',', @contig_status);

  $num_of_false_contigs += scalar(grep {$_ == 0} @contig_status);
  $num_of_contigs += @alleles;

  if(!defined($vars{$vkey})) {
    $vcf->print_entry($entry, $out);
  }
}

print "False contigs: ".pretty_fraction($num_of_false_contigs, $num_of_contigs)."\n";

close($fh);
close($out);

#
# Functions
#

sub any { $_ && return 1 for @_; 0 }
sub all { $_ || return 0 for @_; 1 }

sub contig_exists
{
  my ($contig) = @_;
  $contig =~ s/^N+//i;
  return contig_exists_sub($contig) || contig_exists_sub(rev_comp($contig));
}

sub contig_exists_sub
{
  my ($contig) = @_;
  my ($idx,$pos) = $search_genomes->find_index($contig);
  return ($idx != -1);
  # my $i = 0;
  # while($i < @genomes && index($genomes[$i], $contig) == -1) {$i++;}
  # return ($i < @genomes);
}

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

