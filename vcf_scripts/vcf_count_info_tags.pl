#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

# Config #
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_tag.pl <file.vcf> <infotag1 ..>\n";
  print STDERR "  Prints comma separated histogram of info tag value counts\n";
  print STDERR "  If <file.vcf> is '-' reads from STDIN\n";
  exit;
}

if(@ARGV < 2)
{
  print_usage();
}

my $vcf_file = shift;

my @tags = @ARGV;

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

close($vcf_handle);

# Done.

sub print_set
{
  my ($hashref, @tag_values) = @_;

  if(@tag_values == @tags)
  {
    # print tags
    print join($csvsep, @tag_values) . $csvsep .
          $hashref->{$tag_values[$#tags]} . "\n";
  }
  else
  {
    for my $value (sort keys %{$hashref})
    {
      print_set($hashref, @tag_values, $value);
    }
  }
}