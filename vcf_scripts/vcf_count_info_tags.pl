#!/usr/bin/env perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;

# Config #
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_print_tag.pl <file.vcf> <infotag1 ..>
  Prints comma separated histogram of info tag value counts
  If <file.vcf> is '-' reads from STDIN\n";

  exit(-1);
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 2)
{
  print_usage();
}

my $vcf_file = shift;

my @tags = @ARGV;

#
# Open VCF File
#
my $vcf = vcf_open($vcf_file);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $counts = {};

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my @data = ();

  my $hashref = $counts;

  for(my $i = 0; $i < $#tags; $i++)
  {
    my $value = defined($vcf_entry->{'INFO'}->{$tags[$i]})
                ? $vcf_entry->{'INFO'}->{$tags[$i]} : "";

    if(!defined($hashref->{$value}))
    {
      $hashref->{$value} = {};
    }

    $hashref = $hashref->{$value};
  }

  my $value = defined($vcf_entry->{'INFO'}->{$tags[$#tags]})
              ? $vcf_entry->{'INFO'}->{$tags[$#tags]} : "";

  $hashref->{$value}++;
}

# Now print!
# Print CSV header
print join($csvsep, @tags) . $csvsep . "count\n";
print_set($counts);

$vcf->vcf_close();

# Done.

sub print_set
{
  my ($hashref, @tag_values) = @_;

  if(@tag_values == @tags)
  {
    # print tags
    print join($csvsep, @tag_values) . $csvsep . $hashref . "\n";
  }
  else
  {
    my @keys = keys %{$hashref};
    
    # Check if all numbers
    my $all_nums = 1;
    
    for my $key (@keys)
    {
      if($key !~ /^(?:\d+\.?\d*|\d*\.?\d+)$/)
      {
        $all_nums = 0;
        last;
      }
    }
    
    @keys = $all_nums ? sort {$a <=> $b} @keys : sort @keys;
    
    for my $value (@keys)
    {
      print_set($hashref->{$value}, @tag_values, $value);
    }
  }
}
