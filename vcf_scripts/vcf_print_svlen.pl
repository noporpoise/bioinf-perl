#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

# dev: add merge

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_svlen.pl [--abs|--tag <t>] <max_svlen> [file.vcf]\n";
  print STDERR "  Print SVLENs\n";
  print STDERR "  --tag specify which INFO tag to use\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 5)
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

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $max_svlen = shift;
my $vcf_file = shift;

if($max_svlen !~ /^\d+$/)
{
  print_usage("max_svlen should be a positive integer " .
              "('$max_svlen' not allowed)");
}

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

my $vcf_entry;

my @ins = map {0} 0..$max_svlen;
my @del = map {0} 0..$max_svlen;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $svlen = $vcf_entry->{'INFO'}->{$tag};

  if($abs)
  {
    $svlen = abs($svlen);
  }

  if($svlen >= 0)
  {
    if($svlen <= $max_svlen)
    {
      $ins[$svlen]++;
    }
  }
  else
  {
    # deletion
    $svlen = abs($svlen);
    if($svlen <= $max_svlen)
    {
      $del[$svlen]++;
    }
  }
}

close($vcf_handle);

if(!$abs)
{
  print join("\n", map {"-".$_.",".$del[$_]} reverse(1..$max_svlen))."\n";
}

print join("\n", map {$_.",".$ins[$_]} 0..$max_svlen)."\n";

