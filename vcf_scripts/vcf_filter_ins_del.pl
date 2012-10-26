#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./vcf_filter_ins_del.pl [OPTIONS] [in.vcf]
  Prints variants that passed the filtering.  If [in.vcf] is not passed or is '-'
  then reads from stdin.  Only prints variants with the AA=[0,1] INFO tag

  Options:
   --clean   Only print a deletion where one path is of length 0
   --invert  Print variants that failed
   INS       Only print insertions (uses AA tag)
   DEL       Only print deletions (uses AA tag)
   INDEL     Print polarised insertions AND deletions
   +size     Insertions of a given size
   -size     Deletions of a given size
   abs_size  Indels of a given size\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

my $invert = 0;
my $clean_only = 0;
my $filter_insertions;
my $filter_deletions;
my $filter_size;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-invert?$/i)
  {
    shift;
    $invert = 1;
  }
  elsif($ARGV[0] =~ /^-?-clean$/i)
  {
    shift;
    $clean_only = 1;
  }
  elsif($ARGV[0] =~ /^INS$/i)
  {
    shift;
    $filter_insertions = 1;
  }
  elsif($ARGV[0] =~ /^DEL$/i)
  {
    shift;
    $filter_deletions = 1;
  }
  elsif($ARGV[0] =~ /^INDEL$/i)
  {
    shift;
    $filter_insertions = 1;
    $filter_deletions = 1;
  }
  elsif($ARGV[0] =~ /^([\+\-]?)(\d+)$/i)
  {
    shift;

    if(defined($1))
    {
      $filter_insertions = ($1 eq "+");
      $filter_deletions = ($1 eq "-");
    }
    else
    {
      $filter_insertions = 1;
      $filter_deletions = 1;
    }

    $filter_size = $2;
  
    if($filter_size == 0)
    {
      print_usage("Cannot filter indels of size zero");
    }
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;

if(@ARGV > 0)
{
  print_usage();
}

if(!defined($filter_insertions) && !defined($filter_deletions))
{
  print_usage("Must specify one of: +<size> | -<size> | INS | DEL | INDEL");
}
elsif(!defined($filter_insertions))
{
  $filter_insertions = 0;
}
elsif(!defined($filter_deletions))
{
  $filter_deletions = 0;
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
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

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $num_of_printed = 0;
my $num_of_variants = 0;
my $num_of_missing_aa = 0;

my $vcf_entry;
my $clean_indel;

while(defined($vcf_entry = $vcf->read_entry))
{
  $num_of_variants++;

  my $aa = $vcf_entry->{'INFO'}->{'AA'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};
  my $print = 0;

  if($svlen == 0 ||
     ($clean_only && !defined($clean_indel = vcf_get_clean_indel($vcf_entry))))
  {
    # Do nothing
  }
  elsif(!defined($aa) || ($aa ne "0" && $aa ne "1"))
  {
    $num_of_missing_aa++;
  }
  {
    my $size = ($aa eq "0" ? $svlen : -$svlen);

    $print = ($filter_insertions && $size > 0 ||
              $filter_deletions && $size < 0) &&
             (!defined($filter_size) || abs($svlen) == $filter_size);
  }

  if($print != $invert)
  {
    $num_of_printed++;
    $vcf->print_entry($vcf_entry);
  }
}

my $vars = $clean_only ? "clean " : "";

if($filter_insertions && $filter_deletions)
{
  $vars = "insertions and deletions";
}
elsif($filter_insertions)
{
  $vars = "insertions";
}
else
{
  $vars = "deletions";
}

# Print filtered rate
print STDERR "vcf_filter_ins_del.pl: " .
             pretty_fraction($num_of_printed, $num_of_variants) . " $vars " .
             (defined($filter_size) ? "of size " . $filter_size . "bp " : "") .
             "printed\n";

if($num_of_missing_aa > 0)
{
  print STDERR "vcf_filter_ins_del.pl: " .
               pretty_fraction($num_of_missing_aa, $num_of_variants) . " " .
               "clean indels were missing AA info tag or it was not '0' or '1'\n";
}

close($vcf_handle);
