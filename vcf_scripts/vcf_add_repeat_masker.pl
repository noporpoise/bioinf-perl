#!/usr/bin/perl

use strict;
use warnings;

use VCFFile;
use GeneticsModule;
use IntervalList;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_repeat_masker.pl <rmsk.txt> [in.vcf]\n";
  print STDERR "  Add repeat annoations to variants\n";
  print STDERR "  Adds INFO values to VCF (all comma-separated lists):\n";
  print STDERR "  - rmsk, rmsk_left, rmsk_right\n";
  print STDERR "    => repeat classes variant is in OR " .
               "classes to the left/right of variant\n";
  print STDERR "  - rmsk_left_dist, rmsk_right_dist => distances to the left/right\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $rmsk_file = shift;
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
# Load rmsk.txt file
#

open(RMSK, $rmsk_file) or die("Cannot read rmsk file '$rmsk_file'");

my %num_per_class = ();
my %repeat_elements_by_chr = ();

my $rmsk_line;

while($rmsk_line = <RMSK>)
{
  my @rmsk_cols = split(/\t/, $rmsk_line);

  my $rmsk_chr = get_clean_chr_name($rmsk_cols[5]);
  my $rmsk_start = $rmsk_cols[6];
  my $rmsk_end = $rmsk_cols[7];
  my $rmsk_class = $rmsk_cols[11];

  if($rmsk_class =~ /\?$/)
  {
    next;
  }

  $num_per_class{$rmsk_class}++;

  my $repeat = {'chr' => $rmsk_chr,
                'start' => $rmsk_start,
                'end' => $rmsk_end,
                'class' => $rmsk_class};

  push(@{$repeat_elements_by_chr{$rmsk_chr}}, [$rmsk_start, $rmsk_end+1, $repeat]);
}

close(RMSK);

my %interval_lists = ();

for my $chr (keys %repeat_elements_by_chr)
{
  $interval_lists{$chr} = new IntervalList(@{$repeat_elements_by_chr{$chr}});
}

my @rmsk_repeat_classes = sort {$a cmp $b} keys %num_per_class;

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

my $additional_header = '';

for my $class (@rmsk_repeat_classes)
{
  $additional_header .= '##INFO=<ID=rmsk_'.$class.',Number=0,Type=Flag,' .
                        'Description="In repeat masker element">' . "\n";
}

print vcf_add_to_header($vcf->get_header(), $additional_header);

my $total_num_entries = 0;
my %num_in_repeat_class = ();
my %missing_chrs = ();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'true_POS'};
  my $len = length($vcf_entry->{'true_REF'});

  my $interval_list = $interval_lists{$chr};

  if(!defined($interval_list))
  {
    if(!defined($missing_chrs{$chr}))
    {
      print STDERR "Chromosome '$chr' is missing from rmsk file '$rmsk_file'\n";
      $missing_chrs{$chr} = 1;
    }
    next;
  }

  my ($hits_arr,$left_arr,$right_arr)
    = $interval_list->fetch_nearest($pos, $pos+$len);
  
  # Interval List returns EITHER hits OR nearest to left and right

  if(@$hits_arr > 0)
  {
    my @classes = map {$_->{'class'}} @$hits_arr;
    $vcf_entry->{'INFO'}->{'rmsk'} = join(",",@classes);
  
    for my $class (@classes)
    {
      $num_in_repeat_class{$class}++;
    }
  }
  else
  {
    if(@$left_arr > 0)
    {
      my @classes = map {$_->{'class'}} @$left_arr;
      my $dist = $pos - $left_arr->[0]->{'end'};

      $vcf_entry->{'INFO'}->{'rmsk_left'} = join(",", @classes);
      $vcf_entry->{'INFO'}->{'rmsk_left_dist'} = $dist;
    }

    if(@$right_arr > 0)
    {
      my @classes = map {$_->{'class'}} @$right_arr;
      my $dist = $right_arr->[0]->{'start'} - $pos;

      $vcf_entry->{'INFO'}->{'rmsk_right'} = join(",", @classes);
      $vcf_entry->{'INFO'}->{'rmsk_right_dist'} = $dist;
    }
    
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);

print STDERR "Of $total_num_entries VCF entries:";

for my $class (sort {$a cmp $b} keys %num_in_repeat_class)
{
  print STDERR "  $num_in_repeat_class{$class} are in repeat class $class\n";
}

