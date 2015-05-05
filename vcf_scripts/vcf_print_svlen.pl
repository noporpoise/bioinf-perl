#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use UsefulModule;
use VCFFile;

# config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_svlen.pl [--abs|--tag <t>] [file.vcf]\n";
  print STDERR "  Print size distribution of SVLENs\n";
  print STDERR "  --tag specify which INFO tag to use\n";
  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV > 4)
{
  print_usage();
}

my $abs = 0;
my $tag = 'SVLEN';

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-abs$/i)
  {
    shift;
    $abs = 1;
  }
  elsif($ARGV[0] =~ /^-?-t(ag)?$/i)
  {
    shift;
    $tag = shift;
  }
  else
  {
    last;
  }
}

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

my $num_of_variants = 0;
my $num_of_usable_variants = 0;

my @ins = ();
my @del = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $svlen = $vcf_entry->{'INFO'}->{$tag};

  if(defined($svlen))
  {
    $num_of_usable_variants++;

    if($abs)
    {
      $svlen = abs($svlen);
    }

    if($svlen >= 0)
    {
      $ins[$svlen]++;
    }
    else
    {
      # deletion
      $svlen = abs($svlen);
      $del[$svlen]++;
    }
  }
}

$vcf->vcf_close();

print STDERR "vcf_print_svlen.pl: Printing distribution of tag '$tag'".
             ($abs ? ' (abs value)' : '') . "\n";

print STDERR "vcf_print_svlen.pl: " .
             pretty_fraction($num_of_usable_variants, $num_of_variants)." ".
             "variants had tag '$tag'\n";

# Print header

print "size".$csvsep."count\n";

if(!$abs)
{
  print join("", map { "-".$_.$csvsep.(defined($del[$_]) ? $del[$_] : 0)."\n" }
                   reverse(1..$#del));
}

print join("", map {$_ . $csvsep . (defined($ins[$_]) ? $ins[$_] : 0) . "\n" }
                 1..$#ins);

