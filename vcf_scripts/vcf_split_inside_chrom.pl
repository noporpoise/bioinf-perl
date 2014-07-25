#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;
use RefGenome;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_split_inside_chrom.pl <in.vcf> <out.a.vcf> <out.b.vcf> [ref.files ..]
  Print variants within chroms to out.a.vcf.
  Print variants outside chom to out.b.vcf\n";

  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 4)
{
  print_usage();
}
my $vcf_file = shift;
my $valid_file = shift;
my $invalid_file = shift;
my @ref_files = @ARGV;

if(@ref_files == 0)
{
  push(@ref_files, '-');
}

if(!defined($invalid_file))
{
  print_usage();
}

#
# Read VCF
#
my $vcf = vcf_open($vcf_file);

my ($out_valid, $out_invalid);
open($out_valid, ">$valid_file") or die("Cannot open '$valid_file'");
open($out_invalid, ">$invalid_file") or die("Cannot open '$invalid_file'");

#
# Load reference
#
my $genome = new RefGenome();
$genome->load_from_files(@ref_files);

$vcf->print_header($out_valid);
$vcf->print_header($out_invalid);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $fa_chrom = $genome->guess_chrom_fasta_name($vcf_entry->{'CHROM'});

  if(!defined($fa_chrom) || $vcf_entry->{'POS'} > $genome->get_chr_length($fa_chrom))
  {
    $vcf->print_entry($vcf_entry, $out_invalid);
  }
  else
  {
    $vcf->print_entry($vcf_entry, $out_valid);
  }
}

$vcf->vcf_close();
close($out_valid);
close($out_invalid);
