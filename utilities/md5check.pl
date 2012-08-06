#!/usr/bin/perl

use strict;
use warnings;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./md5check.pl <file1 ..>
  Compares md5sum output\n";
}

if(@ARGV < 1)
{
  print_usage();
}

my %file_hashes = ();

for my $file (@ARGV)
{
  open(FILE, $file) or die("Cannot open file '$file'");

  my $line;
  while(defined($line = <FILE>))
  {
    chomp($line);

    my ($hash, $name) = split("  ", $line);
  
    if(!defined($file_hashes{$name}))
    {
      $file_hashes{$name} = {};
      $file_hashes{$name}->{$hash} = [$file];
    }
    elsif(!defined($file_hashes{$name}->{$hash}))
    {
      $file_hashes{$name}->{$hash} = [$file];
    }
    else
    {
      push(@{$file_hashes{$name}->{$hash}}, $file);
    }
  }

  close(FILE);
}

my @files = sort keys %file_hashes;

for my $file (@files)
{
  my @hashes_for_file = keys %{$file_hashes{$file}};

  if(@hashes_for_file > 1)
  {
    print "Error: " . $file . " [".scalar(@hashes_for_file)."]\n";
    for my $hash (@hashes_for_file)
    {
      print "  $hash [".join(" ", @{$file_hashes{$file}->{$hash}})."]\n";
    }
  }
  else
  {
    # File only has one hash - show the number of md5sum files that list it
    print $hashes_for_file[0] . "  " . $file . " [" .
          scalar(@{$file_hashes{$file}->{$hashes_for_file[0]}}) . "]\n";
  }
}
