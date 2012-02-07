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

  print STDERR
"Usage: ./vcf_add_ancestral.pl <stampy.sam> <min_MAPQ> <outgroup_name> <in.vcf> " .
  "<ancestral_ref.fa1 ..>\n" .
"  Assigns ancestral allele to VCF.  Adds AA and AALEN INFO tags. If <in.vcf> \n" .
"  is '-' reads from stdin\n" .
"  \n" .
"  For variants with both alleles the same length, uses ancestral reference.\n" .
"  For indels, uses mapping of 5'+allele+3' flank to outgroup genome.\n" .
"  * AA = [. -> unknown, 0 -> ref, 1 -> alt allele]\n" .
"  * AALEN = (derived length - ancestral length)\n";
  exit;
}

if(@ARGV < 5)
{
  print_usage();
}

# Ancestral Reference from:
# ftp://ftp.ebi.ac.uk/pub/databases/ensembl/jherrero/ancestral/
#
# It uses the same coordinates as PanTro2
#
# Contains:
# ACGT => confident
# acgt => less confident
# . => unknown
# - => deleted

my $mapping_file = shift;
my $min_mapq = shift;
my $outgroup_name = shift;
my $vcf_file = shift;
my @ancestral_ref_files = @ARGV;

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

my @mapping_columns = qw(name flags chr pos MAPQ CIGAR
                         rnext pnext tlen seq quality);

#
# Load Ancestral Reference
#
print STDERR "vcf_add_ancestral_to_indels.pl: Loading ancestral ref...\n";
my ($anc_ref_hash) = read_all_from_files(@ancestral_ref_files);

my @anc_chrs = keys %$anc_ref_hash;
# looks like: >ANCESTOR_for_chromosome:CHIMP2.1:1:1:229974691:1
for my $key (@anc_chrs)
{
  if($key =~ /ANCESTOR_for_chromosome:CHIMP2.1:(?:chr)?([\dXY]*)([ab]*):/i)
  {
    my $chr = "chr".uc($1).lc($2);
    $anc_ref_hash->{$chr} = $anc_ref_hash->{$key};
    delete($anc_ref_hash->{$key});
    #print STDERR "vcf_add_ancestral.pl: Loaded '$chr'\n";
  }
  else
  {
    #print STDERR "vcf_add_ancestral.pl: Ancestral ref ignored: $key\n";
  }
}

#
# Print VCF header
#
my $tag_description = "Ancestral Allele using outgroup $outgroup_name " .
                      "(.=unknown, 0=reference, 1=alternate allele)";

$vcf->add_header_tag("INFO", "AA", 1, "String", $tag_description);
$vcf->print_header();

# Read VCF entry
my $vcf_entry;

my $num_of_variants = 0;
my $num_of_np = 0;
my $num_of_indels = 0;

my $num_of_np_polarised = 0;
my $num_of_np_ref = 0;

my $num_of_indels_polarised = 0;
my $num_of_indels_ref = 0;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $var_id = $vcf_entry->{'ID'};
  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'true_POS'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};

  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});

  my $aa;

  if($svlen == 0)
  {
    # Nucleotide Polymorphism (NP) - use ancestral allele
    $num_of_np++;

    if(!defined($anc_ref_hash->{$chr}))
    {
      print STDERR "vcf_add_ancestral.pl: Ancestor lacks chromosome '$chr'\n";
      die();
    }

    my $anc_ref = substr($anc_ref_hash->{$chr},
                         $vcf_entry->{'true_POS'}-1, # because perl uses 0-based
                         length($ref_allele));

    my $anc_alt = substr($anc_ref_hash->{$chr},
                         $vcf_entry->{'true_POS'}-1, # because perl uses 0-based
                         length($alt_allele));

    my $ref_matches = ancestral_match($anc_ref, $ref_allele);
    my $alt_matches = ancestral_match($anc_alt, $alt_allele);

    if($ref_matches && !$alt_matches)
    {
      $aa = '0';
      $num_of_np_polarised++;
    }
    elsif(!$ref_matches && $alt_matches)
    {
      $aa = '1';
      $num_of_np_polarised++;
      $num_of_np_ref++;
    }
    else
    {
      $aa = '.';
    }
  }
  else
  {
    # Indel - use stampy mapping
    $num_of_indels++;

    # Read corresponding mapping
    my ($ref_mapping_line, $alt_mapping_line);

    my %ref_mapping;
    my %alt_mapping;

    while(defined($ref_mapping_line = <$mapping_handle>) &&
          defined($alt_mapping_line = <$mapping_handle>))
    {
      chomp($ref_mapping_line);
      chomp($alt_mapping_line);

      if(!defined($ref_mapping_line) || !defined($alt_mapping_line))
      {
        die("Ran out of mappings (var: $var_id)");
      }

      %ref_mapping = ();
      %alt_mapping = ();

      @ref_mapping{@mapping_columns} = split(/\t/, $ref_mapping_line);
      @alt_mapping{@mapping_columns} = split(/\t/, $alt_mapping_line);

      if(lc($ref_mapping{'name'}) eq lc($var_id."_ref"))
      {
        last;
      }
    }

    if(!defined($ref_mapping_line))
    {
      die("vcf_add_ancestral.pl: Couldn't find mapping for var '$var_id'\n");
    }
    elsif(lc($alt_mapping{'name'}) ne lc($var_id."_alt"))
    {
      die("Alt alleles' names don't match " .
          "(var '$var_id' => '$ref_mapping{'name'}')");
    }

    if($ref_mapping{'CIGAR'} eq "*" || $alt_mapping{'CIGAR'} eq "*" ||
       $ref_mapping{'chr'} ne $alt_mapping{'chr'} || 
       abs($ref_mapping{'pos'} - $alt_mapping{'pos'}) > 32)
    {
      # Ref/alt don't map, map to different chromosomes or map too far apart
      # 32 used here because it is the flank size
      $aa = '.';
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

      $aa = get_indel_ancestor($long_cigar, $short_cigar,
                               $long_mapq, $short_mapq,
                               $svlen);
      
      if($aa ne '.')
      {
        $num_of_indels_polarised++;

        if($aa eq '0')
        {
          $num_of_indels_ref++;
        }
      }
    }
  }

  $vcf_entry->{'INFO'}->{'AA'} = $aa;

  if($aa eq '.')
  {
    $vcf_entry->{'INFO'}->{'AALEN'} = '.';
  }
  else
  {
    $vcf_entry->{'INFO'}->{'AALEN'} = ($aa eq '0' ? $svlen : -$svlen);
  }

  $vcf->print_entry($vcf_entry);

  # DEBUG
  #if($num_of_variants > 4) {
  #  die("DEBUG END");
  #}
}

my $percent = sprintf("%.2f", 100 * $num_of_np_polarised / $num_of_np);
print STDERR "vcf_add_ancestral.pl: " .
             num2str($num_of_np_polarised) . "/" . num2str($num_of_np) .
             " ($percent%) NP polarised\n";

$percent = sprintf("%.2f", 100 * $num_of_np_ref / $num_of_np_polarised);
print STDERR "vcf_add_ancestral.pl: of which " .
             num2str($num_of_np_ref) . "/" . num2str($num_of_np_polarised) .
             " ($percent%) were ref\n";

$percent = sprintf("%.2f", 100 * $num_of_indels_polarised / $num_of_indels);
print STDERR "vcf_add_ancestral.pl: " .
             num2str($num_of_indels_polarised) . "/" . num2str($num_of_indels) .
             " ($percent%) indels polarised\n";

$percent = sprintf("%.2f", 100 * $num_of_indels_ref / $num_of_indels_polarised);
print STDERR "vcf_add_ancestral.pl: of which " .
             num2str($num_of_indels_ref) . "/" . num2str($num_of_indels_polarised) .
             " ($percent%) were ref\n";

# Done
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

sub ancestral_match
{
  my ($anc, $allele) = @_;

  my $len = length($anc);
  my $unknowns = 0;

  for(my $i = 0; $i < $len; $i++)
  {
    my $anc_char = substr($anc,$i,1);

    if($anc_char eq ".")
    {
      $unknowns++;
    }
    elsif($anc_char ne substr($allele,$i,1))
    {
      return 0;
    }
  }

  return $unknowns < ($len / 2);
}

