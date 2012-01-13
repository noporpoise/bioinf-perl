#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_filter_slippage.pl [--invert | --flag <flag>] [in.vcf]\n";
  print STDERR "  Prints (or flags) slippage indels\n";
  print STDERR "  Requires INFO tags 'left_flank' and 'right_flank'\n";
  print STDERR "\n";
  print STDERR "  --invert        Print or flag non slippage indels\n";
  print STDERR "  --flag <flag>   Add a given flag\n";
  exit;
}

if(@ARGV > 4)
{
  print_usage();
}

my $invert = 0;
my $flag;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-invert$/i)
  {
    shift;
    $invert = 1;
  }
  elsif($ARGV[0] =~ /^-?-flag$/i)
  {
    shift;
    $flag = shift;
    
    if(!defined($flag))
    {
      print_usage("Flag name needs to be specified");
    }
  }
  elsif(@ARGV > 1)
  {
    print_usage("Unknown argument '$ARGV[0]'");
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
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

if(defined($flag))
{
  my $add_header = "##INFO=<ID=$flag,Number=0,Type=Flag, Description=\"" .
                   ($invert ? "Not slippage indels" : "Slippage indels") .
                   "\">\n";

  print vcf_add_to_header($vcf->get_header(), $add_header);
}
else
{
  print $vcf->get_header();
}

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $ref = $vcf_entry->{'true_REF'};
  my $alt = $vcf_entry->{'true_ALT'};

  my $long_allele = $ref;
  my $short_allele = $alt;

  if(length($long_allele) < length($short_allele))
  {
    # swap
    ($long_allele, $short_allele) = ($short_allele, $long_allele);
  }

  my $is_slippage = 0;

  if(length($short_allele) == 0 && length($long_allele) > 0)
  {
    if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
       !defined($vcf_entry->{'INFO'}->{'right_flank'}))
    {
      print_usage("left_flank, right_flank INFO tags not set");
    }

    my $flanks = substr($vcf_entry->{'INFO'}->{'left_flank'},-length($long_allele)) .
                 substr($vcf_entry->{'INFO'}->{'right_flank'},0,length($long_allele));
    
    $is_slippage = ($flanks =~ /$long_allele/i);
  }

  if(defined($flag))
  {
    # Print all, tag those that are filtered
    if($is_slippage != $invert)
    {
      $vcf_entry->{'INFO_flags'}->{$flag} = 1;
    }

    $vcf->print_entry($vcf_entry);
  }
  elsif($is_slippage != $invert)
  {
    # Just print / filter
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
