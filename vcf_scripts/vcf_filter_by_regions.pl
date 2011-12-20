#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use IntervalList;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_filter_by_regions.pl <OPTIONS> [in.vcf]\n";
  print STDERR "  Filter by regions.  Range is inclusive. " .
               "Does not assume VCF is sorted.  \n";
  print STDERR "  OPTIONS:\n";
  print STDERR "   chr:*          an entire chromosome\n";
  print STDERR "   chr:start-end  a region of a chromosome\n";
  print STDERR "   --file         file containing one arg per line as above\n";
  print STDERR "\n";
  print STDERR "  Examples:\n";
  print STDERR "  \$ ./vcf_regions.pl chr1:10,000-12,000 data.vcf\n";
  print STDERR "  \$ ./vcf_regions.pl --file regions.txt data.vcf\n";
  print STDERR "  \$ cat data.vcf | ./vcf_regions.pl chr1:* chr2:500-1000\n";
  exit;
}

if(@ARGV == 0)
{
  print_usage();
}

my $vcf_file;

my %print_all_chrs = ();
my %regions_by_chr = ();

for(my $i = 0; $i < @ARGV; $i++)
{
  my $arg = $ARGV[$i];

  if(!parse_region($arg))
  {
    if($arg =~ /-?-f(ile)?/i)
    {
      if($i == @ARGV - 1)
      {
        print_usage("Missing file argument after '$arg'");
      }

      my $file = $ARGV[++$i];

      open(FILE, $file) or print_usage("Cannot open region file '$file'");
    
      my $line;
      while(defined($line = <FILE>))
      {
        chomp($line);
        # Check line is not empty or a comment line
        if($line !~ /^\s*$/ && $line !~ /^#/ && !parse_region($line))
        {
          print_usage("Unexpected region in file '$file': '$line'");
        }
      }
    
      close(FILE);
    }
    else
    {
      $vcf_file = $arg;
    
      if($i != @ARGV - 1)
      {
        print_usage("Invalid commandline options\n");
      }
    }
  }
}

# Create intervals
my %interval_lists_by_chr = ();

for my $chr (keys %regions_by_chr)
{
  $interval_lists_by_chr{$chr} = new IntervalList(@{$regions_by_chr{$chr}});
}

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

print $vcf->get_header();

my $num_of_filtered_entries = 0;
my $total_num_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'POS'};
  my $len = length($vcf_entry->{'true_REF'});

  my $print = 0;

  if(defined($print_all_chrs{$chr}))
  {
    $print = 1;
  }
  elsif(defined($interval_lists_by_chr{$chr}))
  {
    my @hits = $interval_lists_by_chr{$chr}->fetch($pos, $pos+$len);
    $print = (@hits > 0);
  }

  if($print)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_filtered_entries++;
  }
}

close($vcf_handle);

print STDERR "$num_of_filtered_entries / $total_num_entries printed\n";

sub parse_region
{
  my ($arg) = @_;

  if($arg =~ /^(.+):(?:([0-9,]+)-([0-9,]+)|\*)$/)
  {
    my $chr = $1;
    my $start = $2;
    my $end = $3;

    if(defined($start))
    {
      # Remove commas
      $start =~ s/,//g;
      $end =~ s/,//g;
    
      if($end < $start)
      {
        print_usage("end position is less than start position in: '$arg'\n");
      }

      if($start <= 0 || $end <= 0)
      {
        print_usage("'$arg' is not a valid region (1-based coords)");
      }

      if(!defined($regions_by_chr{$chr}))
      {
        $regions_by_chr{$chr} = [];
      }

      push(@{$regions_by_chr{$chr}}, [$start, $end+1]);
    }
    else
    {
      $print_all_chrs{$chr} = 1;
    }

    return 1;
  }
  else
  {
    return 0;
  }
}
