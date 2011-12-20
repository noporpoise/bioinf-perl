#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "usage: ./fastn_substr.pl <chr:pos:len,chr:pos-end,.> [in1.fa in2.fq ..]\n";
  print STDERR "  Prints substrings from FASTA/Q files (or STDIN if '-')\n";
  print STDERR "  Takes comma separated list of regions in the form:\n";
  print STDERR "    'chr:start:length,..' and 'chr:start-end,..'\n";
  print STDERR "  Coordinates are 1-based.  " .
               "Start/end positions < 0 mean X bp from the end.\n";
  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $searches_list = shift;
my @files = @ARGV;

my @searches = split(",", $searches_list);
my %search_chrs = ();

for my $search (@searches)
{
  my ($chr,$start,$length,$end);
  
  if($search =~ /^(.*):(-?\d+):(\d+)$/)
  {
    $chr = $1;
    $start = $2;
    $length = $3;
    
    if($start == 0) {
      print_usage("Start position is 1-based - cannot be 0");
    }
    elsif($length == 0) {
      print_usage("Substring length cannot be 0");
    }
    elsif($start < 0 && $length > -$start)
    {
      print_usage("Start position is " . (-$start) . "bp from the end of the " .
                  "read, but length is ".$length."bp (length too long!)");
    }
  }
  elsif($search =~ /^(.*):(-?\d+)-(-?\d+)$/)
  {
    $chr = $1;
    $start = $2;
    $end = $3;
    
    if($start == 0) {
      print_usage("Start position is 1-based - cannot be 0");
    }
    elsif($end == 0) {
      print_usage("End position cannot be 0");
    }
    elsif($start < 0 && $end < $start)
    {
      print_usage("Start position is " . (-$start) . "bp from the end of the " .
                  "read, but end is ".(-$end)."bp from the end (negative length!)");
    }
  }
  else {
    print_usage("Invalid position argument '$search'");
  }

  if(!defined($search_chrs{$chr}))
  {
    # First search position on this chromosome
    $search_chrs{$chr} = [];
  }
  
  push(@{$search_chrs{$chr}}, [$start, $length, $end]);
}


#
# Open FASTA/Q handles
#
if(@files == 0)
{
  my $fastn_handle = open_stdin();
  
  search_file($fastn_handle);
  
  close($fastn_handle);
}
else
{
  for my $file (@files)
  {
    my $fastn_handle;

    if($file eq "-")
    {
      $fastn_handle = open_stdin();
    }
    else
    {
      open($fastn_handle, $file)
        or print_usage("Cannot open FASTA/Q file '$file'");
    }
    
    search_file($fastn_handle);
    
    close($fastn_handle);
    
    if(keys(%search_chrs) == 0)
    {
      # All entries have been found - no need to open any more files
      last;
    }
  }
}

# Done

if(keys(%search_chrs) > 0)
{
  # Print warnings for search results with no match
  for my $chr (sort {$a cmp $b} keys %search_chrs)
  {
    my $array_ref = $search_chrs{$chr};

    for(my $i = 0; $i < @$array_ref; $i++) {
      my $start = $array_ref->[$i]->[0];
      my $length = $array_ref->[$i]->[1];
      print STDERR "# Warning: Couldn't find $chr:$start:$length\n";
    }
  }
}


#
# Functions
#

sub open_stdin
{
  my $stdin_handle;

  if(-p STDIN) {
    # STDIN is connected to a pipe
    open($stdin_handle, "<&=STDIN") or print_usage("Cannot read STDIN pipe");
  }
  else
  {
    print_usage("Must specify or pipe in a FASTA/Q file");
  }
  
  return $stdin_handle;
}

sub search_file
{
  my ($handle) = @_;

  my $fastn = new FASTNFile($handle);

  # Read in
  my ($name, $seq) = $fastn->read_next();

  while(defined($name))
  {
    #print STDERR "Read: '$name' (length ".length($seq).")\n";
  
    my $array_ref = $search_chrs{$name};

    if(defined($array_ref))
    {
      # Print substrings on this chromosome
      for(my $i = 0; $i < @$array_ref; $i++)
      {
        my ($start,$length,$end) = @{$array_ref->[$i]};

        if($start < 0)
        {
          # still in 1-based coords
          $start += length($seq)+1;
        }

        if(!defined($length))
        {
          if($end < 0)
          {
            $end += length($seq)+1;
          }

          $length = $end-$start+1;
          print ">$name:$start-$end\n";
        }
        else
        {
          $end = $start+$length-1;
          print ">$name:$start:$length\n";
        }

        if($end > length($seq)) {
          print STDERR "# Warning: $name:$start:$length is out of bounds of " .
                       "$name:1:".length($seq)."\n";
        }

        # Correct for 1-based coords here (convert to 0-based)
        print substr($seq, $start-1, $length)."\n";
      }
    
      # Only return the first result for each search - so exit
      delete($search_chrs{$name});
    }
  
    if(keys(%search_chrs) == 0)
    {
      last;
    }

    ($name, $seq) = $fastn->read_next();
  }

  close($handle);
}
