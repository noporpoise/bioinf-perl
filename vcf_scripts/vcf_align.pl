#!/usr/bin/perl

use strict;
use warnings;

use GeneticsModule;
use UsefulModule; # num2str
use VCFFile;
use FASTNFile;

use List::Util qw(min max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_align.pl <LEFT|RIGHT> <in.vcf> <ref1.fa ..>\n" .
"  Shift clean indel variants to the left/right.  Variants that do not match\n" .
"  the reference and those that are not clean indels are printed unchanged.\n" .
"  FASTA entry names must match VCF CHROM column.  If <in.vcf> is '-', reads\n".
"  from STDIN.\n";
  exit;
}

if(@ARGV < 3)
{
  print_usage();
}

my $justify = lc(shift);
my $vcf_file = shift;
my @ref_files = @ARGV;

if($justify ne "left" && $justify ne "right")
{
  print_usage("neither 'left' nor 'right' given");
}

#
# Open VCF Handle
#
my $vcf_handle;

if($vcf_file ne "-")
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
# Load reference
#
my ($ref_genomes_hashref, $qual) = read_all_from_files(@ref_files);

my %ref_genomes = %$ref_genomes_hashref;

# Correct '1..' to 'chr1...' etc
# Change chromosome to uppercase
while(my ($key,$value) = each(%ref_genomes))
{
  my $new_key = get_clean_chr_name($key);

  if($new_key ne $key)
  {
    $ref_genomes{$new_key} = uc($ref_genomes{$key});
    delete($ref_genomes{$key});
  }
  else {
    $ref_genomes{$key} = uc($ref_genomes{$key});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Add justify info to header and print
$vcf->add_header_metainfo("variants_justified",$justify);
$vcf->print_header();

my $vcf_entry;

my $num_ref_mismatch = 0;
my $num_of_variants = 0;

my %missing_chrs = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $chr = $vcf_entry->{'CHROM'};
  
  if(!defined($ref_genomes{$chr}))
  {
    $missing_chrs{$chr} = 1;
    next;
  }

  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};

  my $indel = get_clean_indel($vcf_entry);

  my $ref_seq = substr($ref_genomes{$chr}, $var_start, length($ref_allele));

  if($ref_seq ne uc($ref_allele))
  {
    $num_ref_mismatch++;
  }
  elsif(defined($indel))
  {
    if($justify eq "right")
    {
      # justify right
      # while base after variant (on the reference) equals first base of indel
      while(substr($ref_genomes{$chr}, $var_start+length($ref_allele), 1) eq
            substr($indel, 0, 1))
      {
        $indel = substr($indel, 1) . substr($indel, 0, 1);
        $var_start++;
      }
    }
    else
    {
      # justify left
      # while base before variant (on the reference) equals last base of indel
      while(substr($ref_genomes{$chr}, $var_start-1, 1) eq
            substr($indel, -1))
      {
        $indel = substr($indel, -1) . substr($indel, 0, -1);
        $var_start++;
      }
    }

    # Update VCF entry values
    # $var_start was 0-based, VCF POS is 1-based
    $vcf_entry->{'true_POS'} = $var_start+2;
    $vcf_entry->{'POS'} = $var_start+1;
    
    $vcf_entry->{length($ref_allele) > 0 ? 'REF' : 'ALT'} = $indel;

    my $prior_base = substr($ref_genomes{$chr}, $var_start-1, 1);
    $vcf_entry->{'true_REF'} = $prior_base.$vcf_entry->{'REF'};
    $vcf_entry->{'true_ALT'} = $prior_base.$vcf_entry->{'ALT'};
  }

  $vcf->print_entry($vcf_entry);
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_add_flanks.pl: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

if($num_ref_mismatch > 0)
{
  print STDERR "vcf_add_flanks.pl: " .
               pretty_fraction($num_ref_mismatch, $num_of_variants) . " " .
               "variants removed for not matching the reference\n";
}

close($vcf_handle);
