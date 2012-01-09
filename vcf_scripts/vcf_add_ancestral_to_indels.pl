#!/usr/bin/perl

use strict;
use warnings;

use Fcntl qw(SEEK_CUR);

use VCFFile;
use FASTNFile;
use UsefulModule;

use List::Util qw(min max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_ancestral_to_indels.pl " .
               "<stampy.sam> <min_MAPQ> <outgroup_name> [in.vcf]\n";
  print STDERR "  Uses mapping of 5'+allele+3' flank to outgroup genome to " .
               "assign ancestral allele to VCF\n";
  print STDERR "  Adds AA=[.,0,1] INFO tag to <in.vcf>.  " .
               ". = unknown, 0 = ref, 1 = alt allele\n";
  print STDERR "  If <in.vcf> is '-', reads from stdin\n";
  exit;
}

if(@ARGV < 3 || @ARGV > 4)
{
  print_usage();
}

my $mapping_file = shift;
my $min_mapq = shift;
my $outgroup_name = shift;
my $vcf_file = shift;

# Commandline arg checks
if($min_mapq !~ /^\d+$/)
{
  print_usage("Min. MAPQ must be a positive integer (>=0)");
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
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

my $vcf = new VCFFile($vcf_handle);

#
# Open mapping file
#
my $mapping_handle;

open($mapping_handle, $mapping_file)
  or die("Cannot open mapping file '$mapping_file'");

# skip header lines starting @...
my $mapping_line;
while(defined($mapping_line = <$mapping_handle>) && $mapping_line =~ /^@/) {}

seek($mapping_handle,-length($mapping_line),SEEK_CUR);

# Read VCF
my $header_add = "##INFO=<ID=AA,Number=1,Type=String," .
                 "Description=\"Ancestral Allele using outgroup $outgroup_name " .
                 "(.=unknown, 0=reference, 1=alternate allele)\">\n";

print vcf_add_to_header($vcf->get_header(), $header_add);

my @mapping_columns = qw(name flags chr pos MAPQ CIGAR
                         rnext pnext tlen seq quality);

# Read VCF entry
my $vcf_entry;

my $num_of_indels = 0;
my $num_of_variants = 0;
my $num_of_ancestral = 0;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $var_id = $vcf_entry->{'ID'};
  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'true_POS'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};

  if($svlen == 0)
  {
    # Not indel, skip
    # First check if we need to read in any mapping lines
    my $ref_mapping_line = <$mapping_handle>;
    my $alt_mapping_line = <$mapping_handle>;
    
    if(defined($ref_mapping_line) && defined($alt_mapping_line))
    {
      # Remember length in case we need to back track 
      my $read_len = length($ref_mapping_line) + length($alt_mapping_line);

      my %ref_mapping = ();

      @ref_mapping{@mapping_columns} = split(/\t/, $ref_mapping_line);
    
      if(lc($ref_mapping{'name'}) ne lc($var_id."_ref"))
      {
        # backtrack
        seek($mapping_handle, -$read_len, SEEK_CUR);
      }
    }

    next;
  }

  $num_of_indels++;

  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});
  
  my $left_flank = $vcf_entry->{'left_flank'};
  my $right_flank = $vcf_entry->{'right_flank'};

  my $ref_mapping_line = <$mapping_handle>;
  my $alt_mapping_line = <$mapping_handle>;

  chomp($ref_mapping_line);
  chomp($alt_mapping_line);

  if(!defined($ref_mapping_line) || !defined($alt_mapping_line))
  {
    die("Ran out of mappings (var: $var_id)");
  }

  my %ref_mapping = ();
  my %alt_mapping = ();

  @ref_mapping{@mapping_columns} = split(/\t/, $ref_mapping_line);
  @alt_mapping{@mapping_columns} = split(/\t/, $alt_mapping_line);

  # Sanity checks
  if(lc($ref_mapping{'name'}) ne lc($var_id."_ref"))
  {
    die("Reference alleles' names don't match " .
        "(var '$var_id' => '$ref_mapping{'name'}')");
  }
  elsif(lc($alt_mapping{'name'}) ne lc($var_id."_alt"))
  {
    die("Alt alleles' names don't match " .
        "(var '$var_id' => '$ref_mapping{'name'}')");
  }
  
  # 32 used here because it is the flank size
  
  if($ref_mapping{'CIGAR'} eq "*" || $alt_mapping{'CIGAR'} eq "*" ||
     $ref_mapping{'chr'} ne $alt_mapping{'chr'} || 
     abs($ref_mapping{'pos'} - $alt_mapping{'pos'}) > 32)
  {
    $vcf_entry->{'INFO'}->{'AA'} = '.';
  }
  else
  {
    if($ref_mapping{'CIGAR'} !~ /^(?:\d+[MID])+$/i)
    {
      # Full range MIDNSHP=X
      die("Unexpected cigar entry (var: $var_id; cigar: '$ref_mapping{'CIGAR'}')");
    }
    elsif($alt_mapping{'CIGAR'} !~ /^(?:\d+[MID])+$/i)
    {
      die("Unexpected cigar entry (var: $var_id; cigar: '$alt_mapping{'CIGAR'}')");
    }

    my $long_cigar = $ref_mapping{'CIGAR'};
    my $short_cigar = $alt_mapping{'CIGAR'};

    my $long_mapq = $ref_mapping{'MAPQ'};
    my $short_mapq = $alt_mapping{'MAPQ'};

    if($svlen > 0)
    {
      # alt is longer than ref -> swap
      ($short_cigar, $long_cigar) = ($long_cigar, $short_cigar);
      ($short_mapq, $long_mapq) = ($long_mapq, $short_mapq);
    }

    my $aa = get_indel_ancestor($long_cigar, $short_cigar,
                                $long_mapq, $short_mapq,
                                $svlen);

    $vcf_entry->{'INFO'}->{'AA'} = $aa;

    if($aa ne '.')
    {
      $num_of_ancestral++;
    }
  }

  $vcf->print_entry($vcf_entry);
  
  # DEBUG
  #if($num_of_variants > 4) {
  #  die("DEBUG END");
  #}
}

# Check for remaining lines in the mapping file
while(defined($mapping_line = <$mapping_handle>))
{
  chomp($mapping_line);
  if(length($mapping_line) > 0)
  {
    die("Remaining mapping lines: '$mapping_line'");
  }
}

my $percent = 100 * $num_of_ancestral / $num_of_variants;
print STDERR num2str($num_of_ancestral) . " / " . num2str($num_of_variants) .
             " (" . sprintf("%.2f", $percent) . "%) " .
             "polarised using '$outgroup_name'\n";

close($vcf_handle);
close($mapping_handle);


#
# Methods
#

sub get_indel_ancestor
{
  my ($long_cigar, $short_cigar, $longer_mapq, $shorter_mapq, $svlen) = @_;
  my $size = abs($svlen);

  my $long_cigar_match = ($long_cigar =~ /^\d+M$/i);
  my $short_cigar_match = ($short_cigar =~ /^\d+M$/i);

  if($short_cigar_match && !$long_cigar_match && $shorter_mapq >= $min_mapq)
  {
    # shorter allele matches
    # If (alt-ref) > 0 then reference (0) is the shorter allele
    #                  otherwise alt (1)
    return $svlen > 0 ? '0' : '1';
  }
  elsif($long_cigar_match && !$short_cigar_match && $longer_mapq >= $min_mapq)
  {
    # longer allele matches
    # If (alt-ref) > 0 then alt (1) is the longer allele
    #                  otherwise ref (0)
    return $svlen > 0 ? '1' : '0';
  }
  else
  {
    return '.';
  }
}
