#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;

my $ins_tag = 'INSERTION';
my $del_tag = 'DELETION';

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_ins_del_tags.pl [file.vcf]\n";
  print STDERR "  Uses AA tag to label indels in the INFO field\n";
  print STDERR "  Tags with: $ins_tag & $del_tag\n";
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

my $add_header
  = "##INFO=<ID=$ins_tag,Number=1,Type=Flag,Description=\"An insertion\">\n" .
    "##INFO=<ID=$del_tag,Number=1,Type=Flag,Description=\"A deletion\">";

print vcf_add_to_header($vcf->get_header(), $add_header);

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $ins = 0;
  my $set = 0;

  # SVLEN is length(alt) - length(ref)
  if($vcf_entry->{'INFO'}->{'AA'} eq "0")
  {
    # ref allele is ancestral
    $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} > 0);
    $set = 1;
  }
  elsif($vcf_entry->{'INFO'}->{'AA'} eq "1")
  {
    # alt allele is ancestral
    $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} < 0);
    $set = 1;
  }

  if($set)
  {
    if($ins)
    {
      $vcf_entry->{'INFO'}->{$ins_tag} = 1;
      delete($vcf_entry->{'INFO'}->{$del_tag});
    }
    else
    {
      $vcf_entry->{'INFO'}->{$del_tag} = 1;
      delete($vcf_entry->{'INFO'}->{$ins_tag});
    }
  }
  else
  {
    delete($vcf_entry->{'INFO'}->{$del_tag});
    delete($vcf_entry->{'INFO'}->{$ins_tag});
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);
