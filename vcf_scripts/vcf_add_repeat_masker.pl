#!/usr/bin/perl

use strict;
use warnings;

use List::MoreUtils qw(uniq);

use VCFFile;
use GeneticsModule;
use IntervalList;
use UsefulModule; # num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_repeat_masker.pl <rmsk.txt> [in.vcf]\n";
  print STDERR "  Add repeat annoations to variants\n";
  print STDERR "  Adds INFO values to VCF :\n";
  print STDERR "  - rmsk, rmsk_left, rmsk_right (comma-separated lists)\n";
  print STDERR "    repeat classes variant is in OR " . 
               "classes to the left/right of variant\n";
  print STDERR "  - rmsk_left_dist, rmsk_right_dist => distances to the " .
               "left/right\n";
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

  if(!defined($repeat_elements_by_chr{$rmsk_chr}))
  {
    $repeat_elements_by_chr{$rmsk_chr} = {};
  }
  
  if(!defined($repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}))
  {
    $repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class} = [];
  }

  push(@{$repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}},
       [$rmsk_start, $rmsk_end+1, $repeat]);
}

close(RMSK);

my @rmsk_classes = sort {$a cmp $b} keys %num_per_class;

my %interval_lists = ();

for my $rmsk_chr (keys %repeat_elements_by_chr)
{
  for my $rmsk_class (@rmsk_classes)
  {
    $interval_lists{$rmsk_chr}->{$rmsk_class}
      = new IntervalList(@{$repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Add header tags
for my $class (@rmsk_classes)
{
  my $tag_description = "Variant in repeat masker element $class";
  $vcf->add_header_tag("INFO", "rmsk_".$class, 0, "Flag", $tag_description);
}

# Print VCF header
$vcf->print_header();

print vcf_add_to_header($vcf->get_header(), $additional_header);

my $total_num_entries = 0;
my %num_in_repeat_class = ();
my %missing_chrs = ();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};

  my $var_start = $vcf_entry->{'true_POS'};
  my $var_end = $var_start + length($vcf_entry->{'true_REF'});

  my $interval_lists_hashref = $interval_lists{$chr};

  if(!defined($interval_lists_hashref))
  {
    if(!defined($missing_chrs{$chr}))
    {
      print STDERR "Chromosome '$chr' is missing from rmsk file '$rmsk_file'\n";
      $missing_chrs{$chr} = 1;
    }
    next;
  }

  my @hits = ();
  my @hits_left = ();
  my @hits_right = ();

  for my $rmsk_class (@rmsk_classes)
  {
    my ($hits_arr, $left_arr, $right_arr)
      = $interval_lists_hashref->{$rmsk_class}->fetch_nearest($var_start,
                                                              $var_end);

    # Interval List returns EITHER hits OR nearest to left and right
    push(@hits, @$hits_arr);
    push(@hits_left, @$left_arr);
    push(@hits_right, @$right_arr);
  }

  if(@hits > 0)
  {
    my @classes = map {$_->{'class'}} @hits;
    $vcf_entry->{'INFO'}->{'rmsk'} = join(",", @classes);
  
    for my $class (@classes)
    {
      $num_in_repeat_class{$class}++;
    }
  }
  else
  {
    # Unset any existing values
    delete($vcf_entry->{'INFO'}->{'rmsk'});
  }

  if(@hits_left > 0)
  {
    my @classes = map {$_->{'class'}} @hits_left;
    my @dists = map {$var_start - $_->{'end'}} @hits_left;

    $vcf_entry->{'INFO'}->{'rmsk_left'} = join(",", @classes);
    $vcf_entry->{'INFO'}->{'rmsk_left_dist'} = join(",", @dists);
  }
  else
  {
    # Unset any existing values
    delete($vcf_entry->{'INFO'}->{'rmsk_left'});
    delete($vcf_entry->{'INFO'}->{'rmsk_left_dist'});
  }

  if(@hits_right > 0)
  {
    my @classes = map {$_->{'class'}} @hits_right;
    my @dists = map {$_->{'start'} - $var_end} @hits_right;

    $vcf_entry->{'INFO'}->{'rmsk_right'} = join(",", @classes);
    $vcf_entry->{'INFO'}->{'rmsk_right_dist'} = join(",", @dists);
  }
  else
  {
    # Unset any existing values
    delete($vcf_entry->{'INFO'}->{'rmsk_right'});
    delete($vcf_entry->{'INFO'}->{'rmsk_right_dist'});
  }


  $vcf->print_entry($vcf_entry);
}


print STDERR "vcf_add_repeat_masker.pl: of " .
             num2str($total_num_entries) . " VCF entries:\n";

my @sorted_classes = sort {$num_in_repeat_class{$b} <=> $num_in_repeat_class{$a}}
                     keys %num_in_repeat_class;

@sorted_classes = uniq(@sorted_classes, @rmsk_classes);

for my $class (@sorted_classes)
{
  my $count = $num_in_repeat_class{$class};

  if(!defined($count))
  {
    $count = 0;
  }
  
  my $percent = 100 * $count / $total_num_entries;

  print STDERR "  " . num2str($count) .
               " (" . sprintf("%.2f", $percent) . "%) " .
               " are in repeat class $class\n";
}

close($vcf_handle);
