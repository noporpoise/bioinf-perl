#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max reduce);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # for pretty_fraction

use constant {DUPE_SELECT_NONE => 0,
              DUPE_SELECT_LOWEST_TAG => 1,
              DUPE_SELECT_HIGHEST_TAG => 2,
              DUPE_SELECT_FIRST => 3,
              DUPE_SELECT_LAST => 4};

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_remove_dupes.pl [OPTIONS] [file.vcf]
  Remove entries that match position, REF and ALT alleles. Assumes sorted VCF.

  Order duplicates by an INFO tag and take the variant with highest/lowest value:
  --take_lowest <tag>  OR
  --take_highest <tag> OR
  --take_first         OR
  --take_last          OR
  
  --filter_text <txt>  Add to / set the filter column instead of removing\n";

  exit;
}

if(@ARGV > 5)
{
  print_usage();
}

my $select_dupe;
my $select_tag;

my $filter_txt;

while(@ARGV >= 1)
{
  if($ARGV[0] =~ /^--take_lowest$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_LOWEST_TAG;
    $select_tag = shift
      or print_usage("Missing INFO tag argument to --take_lowest");
  }
  elsif($ARGV[0] =~ /^--take_highest$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_HIGHEST_TAG;
    $select_tag = shift
      or print_usage("Missing INFO tag argument to --take_highest");
  }
  elsif($ARGV[0] =~ /^--take_first$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_FIRST;
  }
  elsif($ARGV[0] =~ /^--take_last$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_LAST;
  }
  elsif($ARGV[0] =~ /^--filter_txt$/i)
  {
    shift;
    $filter_txt = shift or print_usage("Missing FILTER text for --filter_txt");
  }
  elsif(@ARGV > 1)
  {
    print_usage("Unknown option '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'");
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
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

$vcf->print_header();

my $num_of_entries = 0;
my $num_of_printed = 0;

my @prev_variants = ();

my $curr_variant = $vcf->read_entry();

if(defined($curr_variant))
{
  $num_of_entries++;
  push(@prev_variants, $curr_variant);
}

while(defined($curr_variant = $vcf->read_entry()))
{
  $num_of_entries++;

  if($curr_variant->{'CHROM'} ne $prev_variants[0]->{'CHROM'} ||
     $curr_variant->{'POS'} != $prev_variants[0]->{'POS'})
  {
    print_prev_variants();
    @prev_variants = ();
  }

  push(@prev_variants, $curr_variant);
}

print_prev_variants();

print STDERR "vcf_remove_dupes.pl: " .
      pretty_fraction($num_of_printed, $num_of_entries) . " variants printed\n";

close($vcf_handle);



sub print_prev_variants
{
  if(@prev_variants == 1)
  {
    print_variants(@prev_variants);
  }
  else
  {
    vcf_sort_variants(\@prev_variants);

    my @dupes = ($prev_variants[0]);

    for(my $i = 1; $i < @prev_variants; $i++)
    {
      if(uc($prev_variants[$i]->{'ALT'}) ne uc($dupes[0]->{'ALT'}))
      {
        print_variants(@dupes);
        @dupes = ();
      }

      push(@dupes, $prev_variants[$i]);
    }

    print_variants(@dupes);
  }
}

sub print_variants
{
  if(@_ == 1)
  {
    # Not dupe
    $vcf->print_entry($_[0]);
    $num_of_printed++;
  }
  elsif(defined($select_dupe))
  {
    # Select the variant with the highest/lowest INFO field value
    my $selected_variant;

    if($select_dupe == DUPE_SELECT_HIGHEST_TAG)
    {
      $selected_variant = reduce { $a->{'INFO'}->{$select_tag} >=
                                   $b->{'INFO'}->{$select_tag} ? $a : $b } @_;
    }
    elsif($select_dupe == DUPE_SELECT_LOWEST_TAG)
    {
      $selected_variant = reduce { $a->{'INFO'}->{$select_tag} <=
                                   $b->{'INFO'}->{$select_tag} ? $a : $b } @_;
    }
    elsif($select_dupe == DUPE_SELECT_FIRST)
    {
      $selected_variant = $_[0];
    }
    elsif($select_dupe == DUPE_SELECT_LAST)
    {
      $selected_variant = $_[$#_];
    }

    if(defined($filter_txt))
    {
      # Print all variants, labelling all but the selected one in FILTER column
      for my $variant (@_)
      {
        if($selected_variant != $variant)
        {
          vcf_add_filter_txt($variant, $filter_txt);
        }

        $vcf->print_entry($variant);
        $num_of_printed++;
      }
    }
    else
    {
      # Only print selected variant
      $vcf->print_entry($selected_variant);
      $num_of_printed++;
    }
  }
  elsif(defined($filter_txt))
  {
    # Print all variants, but labelled in the FILTER column
    for my $variant (@_)
    {
      vcf_add_filter_txt($variant, $filter_txt);
      $vcf->print_entry($variant);
      $num_of_printed++;
    }
  }
  else
  {
    # Print only the first variant
    $vcf->print_entry($_[0]);
    $num_of_printed++;
  }
}
