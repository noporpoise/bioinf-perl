#!/usr/bin/env perl

use strict;
use warnings;

sub print_usage
{
  if(@_) {
    print STDERR map {"Error: $_\n"} @_;
  }
  print STDERR "usage: ./vcf_sample_names.pl <names.txt> [input.vcf]\n";
  print STDERR "  where names.txt is a file with two tab separated columns:\n";
  print STDERR "    old_name new_name\n";
  print STDERR "  Prints new header to STDOUT\n";
  exit(-1);
}

if(@ARGV < 1 || @ARGV > 2) { print_usage(); }

my ($names_file, $vcf_file) = @ARGV;
my ($line, $vcf_handle);

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else {
  print_usage("Must specify or pipe in a VCF file");
}

my $header = "";

while(defined($line = <$vcf_handle>) && $line =~ /^##/) {
  $header .= $line;
}

close($vcf_handle);

my $names = `cat $names_file` or die("Cannot read $ARGV[1]");

# Regex '/m' for multiline is needed so ^ and $ match line start/end
my ($hline) = ($header =~ /^(#CHROM.*)$/m);

while($names =~ /^(\w*)\s(\w*)/gm) {
  my ($oldname,$newname) = ($1,$2);
  $hline =~ s/\s$oldname(\b)/\t$newname$1/g;
}

$header =~ s/^#CHROM.*/$hline/m;
print $header;
