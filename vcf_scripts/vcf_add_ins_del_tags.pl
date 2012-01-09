#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_ins_del_tags.pl [file.vcf]\n";
  print STDERR "  Adds INDEL=INS|DEL INFO tag using AA and SVLEN tags\n";
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

my $add_header = "##INFO=<ID=INDEL,Number=1,Type=String,Description=\"" .
                 "INS=insertion, DEL=deletion or .=unknown\">";

print vcf_add_to_header($vcf->get_header(), $add_header);

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # SVLEN is length(alt) - length(ref)
  if($vcf_entry->{'INFO'}->{'AA'} eq "0")
  {
    # ref allele is ancestral
    my $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} > 0);
    $vcf_entry->{'INFO'}->{'INDEL'} = $ins ? 'INS' : 'DEL';
  }
  elsif($vcf_entry->{'INFO'}->{'AA'} eq "1")
  {
    # alt allele is ancestral
    my $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} < 0);
    $vcf_entry->{'INFO'}->{'INDEL'} = $ins ? 'INS' : 'DEL';
  }
  else
  {
    $vcf_entry->{'INFO'}->{'INDEL'} = '.';
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);