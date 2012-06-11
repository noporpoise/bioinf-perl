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
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_correct_strand.pl [--remove_mismatches] <in.vcf> <ref1.fa ..>
  If the REF allele doesn't match the reference but the ALT allele does, switch
  them and correct the genotypes.
  
  --remove_mismatches  Remove mismatches that cannot be fixed\n";

  exit;
}

if(@ARGV < 2)
{
  print_usage();
}

my $remove_ref_mismatch = 0;

if($ARGV[0] =~ /^-?-remove_mismatches$/i)
{
  shift;
  $remove_ref_mismatch = 1;
}

my $vcf_file = shift;
my @ref_files = @ARGV;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN)
{
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
my $genome = new RefGenome(uppercase => 1);
$genome->load_from_files(@ref_files);

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

$vcf->print_header();

my @samples = $vcf->get_list_of_sample_names();

my $vcf_entry;

my $num_of_variants = 0;
my $num_of_ref_mismatch = 0;
my $num_of_switched = 0;

my %missing_chrs = ();

# When the reference allele is a deletion (0bp long allele) there is no way it cannot match the genome

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $chr = $vcf_entry->{'CHROM'};

  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;

  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});

  my $max_allele_len = max(length($ref_allele), length($alt_allele));

  my $genome_seq = $genome->get_chr_substr($chr, $var_start, $max_allele_len);

  my $unfixed_mismatch = 0;

  if(!defined($genome_seq))
  {
    $missing_chrs{$chr} = 1;
  }
  else
  {
    my $genome_seq_ref = substr($genome_seq, 0, length($ref_allele));
    my $genome_seq_alt = substr($genome_seq, 0, length($alt_allele));

    if($ref_allele ne $genome_seq_ref)
    {
      # Mismatch
      $num_of_ref_mismatch++;

      if($alt_allele eq $genome_seq_alt)
      {
        # Switch!
        $num_of_switched++;

        # Switch alleles
        my $tmp = $vcf_entry->{'REF'};
        $vcf_entry->{'REF'} = $vcf_entry->{'ALT'};
        $vcf_entry->{'ALT'} = $tmp;

        # Switch genotypes
        for my $sample (@samples)
        {
          my $genotype = $vcf_entry->{$sample};
          
          if($genotype =~ /^([01])([\/\|])([01])(.*)$/)
          {
            my $gt1 = $1 == 0 ? 1 : 0;
            my $gt2 = $3 == 0 ? 1 : 0;

            $vcf_entry->{$sample} = $gt1.$2.$gt2.$4;
          }
          else
          {
            die("Cannot parse genotype '$genotype' (id: $vcf_entry->{'ID'})");
          }
        }
      }
      else
      {
        # Couldn't fix mismatch
        $unfixed_mismatch = 1;
      }
    }
  }

  if(!$unfixed_mismatch || $remove_ref_mismatch)
  {
    $vcf->print_entry($vcf_entry);
  }
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_align.pl Warning: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

print STDERR "vcf_align.pl: " .
             pretty_fraction($num_of_ref_mismatch, $num_of_variants) . " " .
             "variants did not match the reference " .
             ($remove_ref_mismatch ? "(removed)" : "(left in)") . "\n";

print STDERR "vcf_align.pl: " .
             pretty_fraction($num_of_switched,
                             $num_of_ref_mismatch) . " " .
             "mismatches were switched successfully\n";

close($vcf_handle);
