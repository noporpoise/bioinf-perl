#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use BLASTFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./blast_filter.pl <in.blast> <num_matches> [+require1 ..] [-skip1 ...]
  Filter out entries with match in the top <num_matches>\n";
  exit(-1);
}

if(@ARGV < 3) { print_usage(); }

my $blast_file = shift;
my $num_entries = shift;
my @skips = ();
my @reqs = ();

for my $arg (@ARGV)
{
  my $first = substr($arg, 0, 1);
  my $rest = substr($arg, 1);
  if($first eq "+") { push(@reqs, $rest); }
  elsif($first eq "-") { push(@skips, $rest); }
  else {print_usage("Unknown arg: $arg");}
}

print STDERR "reqs: @reqs\n";
print STDERR "skip: @skips\n";

#
# Open BLAST Handle
#
my $blast_handle;

if($blast_file ne "-") {
  open($blast_handle, $blast_file) or die("Cannot open BLAST file '$blast_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($blast_handle, "<&=STDIN") or die("Cannot read pipe");
}
else {
  print_usage("Must specify or pipe in a BLAST file");
}

#
# Read BLAST
#
my $blast = new BLASTFile($blast_handle);

my @entries;

while((@entries = $blast->read_blast_entries()) &&
      @entries > 0 && defined($entries[0]))
{
  my $count = min($num_entries, @entries-1);
  my $do_print = (@reqs == 0);

  # Check for criteria for printing
  for(my $i = 1; $i <= $count && !$do_print; $i++) {
    for my $req (@reqs) {
      if($entries[$i]->{'chr'} =~ /$req/i) {
        $do_print = 1;
        last;
      }
    }
  }

  # Now filter out using skips
  for(my $i = 1; $i <= $count && $do_print; $i++) {
    for my $skip (@skips) {
      if($entries[$i]->{'chr'} =~ /$skip/i) {
        $do_print = 0;
        last;
      }
    }
  }

  if($do_print) { print $entries[0]; }
}

close($blast_handle);
