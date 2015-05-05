#!/usr/bin/env perl

use strict;
use warnings;

use List::Util qw(max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use RefGenome;
use VCFFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_filter_by_ref.pl [OPTIONS] <in.vcf> <ref1.fa ..>
  Filter out variants that do not match the reference.  If <in.vcf> is '-' then
  read from STDIN.

  OPTIONS:
    --tag <t>       Set the FILTER column of mismatches
    --filter_n <k>  Filter out variants with a non-ACGT base within <k> bp\n";

  exit(-1);
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 1)
{
  print_usage();
}

my $filter_out_ns = 0;
my $filter_dist;
my $filter_tag;

while(@ARGV > 2)
{
  if($ARGV[0] =~ /^-?-filter_n$/i)
  {
    shift;
    $filter_dist = shift;

    if(!defined($filter_dist) || $filter_dist !~ /^\d+$/)
    {
      print_usage("--filter_n <k> must be given a positive integer");
    }
  }
  elsif($ARGV[0] =~ /^-?-tag$/i)
  {
    shift;
    if(!defined($filter_tag = shift))
    {
      print_usage("--tag <txt>  requires a tag");
    }
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;
my @ref_files = @ARGV;

if(@ref_files == 0)
{
  print_usage("No references files given");
}


#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

#
# Load reference files
#
my $genome = new RefGenome(uppercase => 1);
$genome->load_from_files(@ref_files);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

# Things read in
my $total_num_entries = 0;

# Things filtered out
my $num_of_ref_mismatches = 0;
my $num_of_filtered_regions = 0;

# Thing printed out
my $num_of_filtered_entries = 0;
my $num_of_printed_entries = 0;

my $vcf_entry;
my %chrom_not_in_ref = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $print = 1;

  my $chr = $vcf_entry->{'CHROM'};

  if(!$genome->chr_exists($chr))
  {
    $print = 0;

    if(!defined($chrom_not_in_ref{$chr}))
    {
      $chrom_not_in_ref{$chr} = 1;

      print STDERR "vcf_filter_by_ref.pl - Warning: chromosome '$chr' not in " .
                   "reference, filtering out variants\n";
    }
  }
  else
  {
    my $ref_allele = uc($vcf_entry->{'REF'});

    # $var_start is the 0-based position of the padding base(s)
    my $var_start = $vcf_entry->{'POS'} - 1;
    my $var_length = length($ref_allele);

    my $chrom_length = $genome->get_chr_length($chr);

    if($var_start + $var_length > $chrom_length)
    {
      $print = 0;
      print STDERR "vcf_filter_by_ref.pl - Warning: variant " .
                   $vcf_entry->{'ID'} . " " .
                   "[".$chr.":".$vcf_entry->{'POS'}.":".$var_length."] " .
                   "out of bounds of ref " .
                   "'".$genome->guess_chrom_fasta_name($chr)."' " .
                   "[length:$chrom_length]\n";
    }
    else
    {
      my $ref_seq = $genome->get_chr_substr0($chr, $var_start, $var_length);

      if($ref_seq ne $ref_allele)
      {
        $print = 0;
        $num_of_ref_mismatches++;
      }
      elsif($filter_out_ns)
      {
        my $region_start = max(0, $var_start - $filter_dist + 1);
        my $region_end = min($chrom_length, $var_start+$var_length+$filter_dist);
        my $region_len = $region_end - $region_start;

        my $region = $genome->get_chr_substr0($chr, $region_start, $region_len);

        if($region =~ /[^acgt]/i)
        {
          $print = 0;
          $num_of_filtered_regions++;
        }
        else
        {
          $num_of_filtered_entries++;
        }
      }
      else
      {
        $num_of_filtered_entries++;
      }
    }
  }

  if($print || defined($filter_tag))
  {
    if(!$print)
    {
      vcf_add_filter_txt($vcf_entry, $filter_tag);
    }

    $num_of_printed_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

print STDERR "vcf_filter_by_ref.pl: " .
             pretty_fraction($num_of_ref_mismatches, $total_num_entries) .
             " variants mismatched\n";

if($filter_out_ns)
{
  print STDERR "vcf_filter_by_ref.pl: " .
               pretty_fraction($num_of_filtered_regions, $total_num_entries) .
               " variants had non-bases within regions of $filter_dist bp\n";
}

print STDERR "vcf_filter_by_ref.pl: " .
             pretty_fraction($num_of_filtered_entries, $total_num_entries) .
             " variants passed filters\n";

print STDERR "vcf_filter_by_ref.pl: " .
             pretty_fraction($num_of_printed_entries, $total_num_entries) .
             " variants printed\n";

$vcf->vcf_close();
