#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use UsefulModule; # for num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./vcf_get_clean_indels.pl [in.vcf]\n";
  print STDERR "  Prints variants were one allele is 0bp and the other is >0bp\n";
  exit;
}

if(@ARGV > 1 || @ARGV == 1 && $ARGV[0] =~ /-?-h(elp)?/i) {
  print_usage();
}

my $vcf_file = shift;

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

my $add_header = "##INFO=<ID=AA_INDEL,Number=1,Type=String," .
                 "Description=\"Whether indel is ins/del using outgroup\">\n";
$add_header .= "##ALT=<ID=INSERT,Description=\"Insertion\">\n";
$add_header .= "##ALT=<ID=DELETE,Description=\"Deletion\">\n";

print vcf_add_to_header($vcf->get_header(), $add_header);

my $num_of_variants = 0;
my $num_of_clean = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  if(length($vcf_entry->{'true_REF'}) == 0 ||
     length($vcf_entry->{'true_ALT'}) == 0)
  {
    my $ancestor = $vcf_entry->{'INFO'}->{'AA'};
    # $ancestor is 0 for ref allele, 1 for alt allele

    if(defined($ancestor) && $ancestor ne ".")
    {
      # Add clearer ancestral annotation to variants
      # (only meaningful on clean indels)

      my $svlen = length($vcf_entry->{'true_ALT'}) -
                  length($vcf_entry->{'true_REF'});

      if($ancestor == 0 && $svlen > 0 ||
         $ancestor == 1 && $svlen < 0)
      {
        $vcf_entry->{'INFO'}->{'AA_INDEL'} = 'INSERT';
      }
      
      if($ancestor == 0 && $svlen < 0 ||
         $ancestor == 1 && $svlen > 0)
      {
        $vcf_entry->{'INFO'}->{'AA_INDEL'} = 'DELETE';
      }
    }

    $num_of_clean++;

    $vcf->print_entry($vcf_entry);
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_of_clean / $num_of_variants;

print STDERR "vcf_get_clean_indels.pl: " .
             num2str($num_of_clean) . " / " . num2str($num_of_variants) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
