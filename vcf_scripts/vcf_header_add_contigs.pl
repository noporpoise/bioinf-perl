#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(first);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;
use lib $FindBin::Bin . '/../lib';

use VCFFile;
use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_header_add_contigs.pl [--assembly <A>] <in.vcf> [<in1.fa> ..]
  Add contig tags to VCF header using FASTA ref files e.g.
    ##contig=<ID=chr1,length=1000000,assembly=None>\n";

  exit(-1);
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV == 0) { print_usage(); }

my $assembly = "None";

if($ARGV[0] =~ /^-?-assembly$/i) { shift; $assembly = shift; }
if(!defined($assembly) || @ARGV < 1) { print_usage(); }

my $vcf_file = shift;
my @ref_paths = @ARGV;

my $vcf = vcf_open($vcf_file);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

#
# Load reference contigs
#
my @contigs = ();
for my $ref_path (@ref_paths) {
  my $fastn = open_fastn_file($ref_path);
  my ($title,$seq);
  while((($title,$seq) = $fastn->read_next()) && defined($title)) {
    $title =~ s/^\s*(\S+).*$/$1/gi;
    push(@contigs, "<ID=$title,length=".length($seq).",Assembly=$assembly>");
  }
  $fastn->close_fastn_file();
}

for my $contig (@contigs) {
  $vcf->add_header_metainfo("contig", $contig);
}

$vcf->print_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $vcf->print_entry($vcf_entry);
}

$vcf->vcf_close();
