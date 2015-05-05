#!/usr/bin/env perl

use strict;
use warnings;

use File::Path qw(make_path remove_tree);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_sort.pl <order> [in.vcf]\n";
  print STDERR "  Order is a comma-separated list of chroms\n";
  print STDERR "  Prints output to STDOUT, prints stats to STDERR\n";
  exit(-1);
}

if(@ARGV < 1 || @ARGV > 2) { print_usage(); }

my $order_list = shift;
my $vcf_file = shift;
my $vcf_handle;

if(!defined($vcf_file))
{
  # open STDIN if it is connected to a pipe
  if(-p STDIN) {
   open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
  }
  else { print_usage("Must specify or pipe in a VCF file"); }
} else {
  open($vcf_handle, $vcf_file) or die("Cannot read VCF: $vcf_file");
}

my @order = split(',', $order_list);

if(grep {$_ =~ /:/} @order) { print_usage("Colons not allowed in chr names");}
if(grep {$_ eq ""} @order) { print_usage("Empty chromosome names not allowed");}

# check for duplicates in order
my %order_hsh = ();
@order_hsh{@order} = 1;
if(scalar(keys(%order_hsh)) < @order) { print_usage("Duplicates in <order>"); }
%order_hsh = ();

# make tmp directory
my $tmpdir = "vsort";
my $rnd;

for(my $i = 0; $i < 1000; $i++ ) {
  $rnd = int(rand(1000000));
  if(!(-e $tmpdir.$rnd)) { last; }
}

$tmpdir .= $rnd;
make_path($tmpdir) or die("Couldn't create tmpdir $tmpdir");

# open chr files in tmp dir
for my $chr (@order) {
  my $file = "$tmpdir/$chr";
  open($order_hsh{$chr}, ">$file") or die("Couldn't open $file");
}

# read header
my $header = "";
my $line;

while(defined($line = <$vcf_handle>) && $line =~ /^#/) {
  $header .= $line;
}

# read and print lines into chrom tmp files
while(defined($line)) {
  my ($chr) = ($line =~ /^(.*?)\t/);
  if(!defined($chr)) { die("VCF format error: $line"); }
  if(defined($order_hsh{$chr})) { my $fh = $order_hsh{$chr}; print $fh $line; }
  $line = <$vcf_handle>;
}

close($vcf_handle);

# close chr files
for my $fh (values %order_hsh) { close($fh); }

# sort separate chromosomes
# merge and print
my %stats = ();
print $header;
for my $chr (@order) {
  my $cmd = "sort -k2 -n $tmpdir/$chr";
  my $num = 0;
  open(CMD, "-|", "$cmd") or die("sort failed: $cmd");
  while (<CMD>) { $num++; print; }
  close(CMD);
  print STDERR "$chr\t$num\n";
}

# delete tmp directory
remove_tree($tmpdir);
