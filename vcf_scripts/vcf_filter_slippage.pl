#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_flank_overlap.pl [in.vcf]\n";
  print STDERR "  requires INFO tags left_flank and right_flank (must be same length)\n";
  print STDERR "  adds INFO tag FO\n";
  exit;
}

if(@ARGV > 1)
{
  print_usage();
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

my $add_header = "##INFO=<ID=FO,Number=1,Type=Integer," .
                 "Description=\"Max number of overlapping bp between flanks\">\n";

print vcf_add_to_header($vcf->get_header(), $add_header);

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $vcf_entry->{'INFO'}->{'FO'}
    = get_flank_overlap($vcf_entry->{'INFO'}->{'left_flank'},
                        $vcf_entry->{'INFO'}->{'right_flank'});

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);

sub get_flank_overlap
{
  my ($left_flank, $right_flank) = @_;

  my $len = length($left_flank);

  for(my $i = 0; $i < $len; $i++)
  {
    if(substr($left_flank,$i) eq substr($right_flank,0,$len-$i))
    {
      return $i;
    }
  }
}
