#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule; # num2str
use VCFFile;
use RefGenome;

# Config
my $multiallelic_tag = 'MULTIALLELIC';
my $dupallele_tag = 'DUP_ALLELE';
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_combine_alleles.pl [in.vcf]
  Combines variants (and their alleles) that have matching ref sites.
  Prints to stdout.  If no in.vcf given, or '-', reads from STDIN\n";

  exit;
}

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
elsif(-p STDIN)
{
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

$vcf->add_header_tag("FILTER", $multiallelic_tag, 0, undef,
                     "Variant formed from merged alleles");

$vcf->add_header_tag("FILTER", $dupallele_tag, 0, undef,
                     "Allele has been merged");

$vcf->print_header();

my @samples = $vcf->get_list_of_sample_names();

my $vcf_entry;

my @vars_at_same_pos = ();

# Stats
my $num_of_entries = 0;
my $num_of_printed = 0;
my $num_of_merges = 0;
my $num_of_vars_merged = 0;
my $biggest_merge = 0;

# Read first variant
if(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_entries++;
  @vars_at_same_pos = ($vcf_entry);
}

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  if($vcf_entry->{'CHROM'} ne $vars_at_same_pos[0]->{'CHROM'} ||
     $vcf_entry->{'POS'} != $vars_at_same_pos[0]->{'POS'})
  {
    # Parse current list
    print_list();

    @vars_at_same_pos = ();
  }

  push(@vars_at_same_pos, $vcf_entry);
}

# Parse remaining list
print_list();

# Print stats

my $mean_vars_per_merge
  = ($num_of_merges == 0 ? 0 : ($num_of_vars_merged / $num_of_merges));

print STDERR "vcf_combine_alleles.pl: mean number of variants per merge is " .
             num2str($mean_vars_per_merge) . "\n";

print STDERR "vcf_combine_alleles.pl: biggest merge has " .
             num2str($biggest_merge) . " variants\n";

print STDERR "vcf_combine_alleles.pl: " .
             pretty_fraction($num_of_printed, $num_of_entries) .
             " variants printed\n";

print STDERR "vcf_combine_alleles.pl: " .
             pretty_fraction($num_of_merges, $num_of_printed) .
             " printed variants are merges\n";

close($vcf_handle);





# Merge and print list of variants starting at the same position
sub print_list
{
  if(@vars_at_same_pos == 1)
  {
    $vcf->print_entry($vars_at_same_pos[0]);
    $num_of_printed++;
  }
  else
  {
    # Sort by length of ref allele
    @vars_at_same_pos = sort {length($a->{'REF'}) <=> length($b->{'REF'})}
                             @vars_at_same_pos;

    for(my $i = 0; $i < @vars_at_same_pos; $i++)
    {
      if($i+1 == @vars_at_same_pos ||
         length($vars_at_same_pos[$i]->{'REF'}) !=
           length($vars_at_same_pos[$i+1]->{'REF'}))
      {
        $vcf->print_entry($vars_at_same_pos[$i]);
        $num_of_printed++;
      }
      else
      {
        # Two variants [$i,$i+1] share the same ref allele - look for more
        my $start = $i;
        my $end = $i+1;

        $num_of_merges++;

        while($end+1 < @vars_at_same_pos &&
              length($vars_at_same_pos[$end]->{'REF'}) ==
                length($vars_at_same_pos[$end+1]->{'REF'}))
        {
          $end++;
        }

        # @vars_at_same_pos[$start..$end] all share the same ref site
        # merge and print a single entry
        my @alleles_arr = map {$_->{'ALT'}} @vars_at_same_pos[$start..$end];
        my %alleles_hash = ();
        @alleles_hash{@alleles_arr} = 1;
        @alleles_arr = sort keys %alleles_hash;

        # May only be one alt allele (duplicated variant)...
        # in which case just print the variants

        if(@alleles_arr == 1)
        {
          for(my $var = $start; $var <= $end; $var++)
          {
            $vcf->print_entry($vars_at_same_pos[$var]);
            $num_of_printed++;
          }
        }
        else
        {
          # Store allele numbers in this hash
          @alleles_hash{$vars_at_same_pos[$end]->{'REF'}} = 0;

          for(my $indx = 1; $indx <= @alleles_arr; $indx++)
          {
            @alleles_hash{$alleles_arr[$indx-1]} = $indx;
          }

          my $alts = join(",", @alleles_arr);

          # Create new variant entry for merge
          my %new_entry_cpy = %{$vars_at_same_pos[$start]};
          my $new_entry = \%new_entry_cpy;

          $new_entry->{'ID'} .= "_merge";
          $new_entry->{'ALT'} = $alts;

          for my $sample (@samples)
          {
            $new_entry->{$sample} = {};
          }

          # Get ploidy
          my $ploidy;
          for(my $var = $start; !defined($ploidy) && $var <= $end; $var++)
          {
            $ploidy = $vcf->get_ploidy($vars_at_same_pos[$var]);
          }

          if(!defined($ploidy))
          {
            $ploidy = 2;
          }

          $new_entry->{'FORMAT'} = ['GT'];

          # For each sample, pick GT with highest conf
          for my $sample (@samples)
          {
            my $max_gt_conf_var;

            for(my $var = $start; $var <= $end; $var++)
            {
              my $gt_conf = $vars_at_same_pos[$var]->{$sample}->{'GT_CONF'};
              my $sample_gt = $vars_at_same_pos[$var]->{$sample}->{'GT_CONF'};

              if(defined($gt_conf) && defined($sample_gt) && $gt_conf ne ".")
              {
                if(!defined($max_gt_conf_var) ||
                   $max_gt_conf_var->{$sample}->{'GT_CONF'}
                     > $vars_at_same_pos[$var]->{$sample}->{'GT_CONF'})
                {
                  $max_gt_conf_var = $vars_at_same_pos[$var];
                }
              }
            }

            if(defined($max_gt_conf_var))
            {
              my $new_gt_call = "";
              my $prev_gt = $max_gt_conf_var->{$sample}->{'GT'};

              # Get alts
              my @alts = split(",", $max_gt_conf_var->{'ALT'});
              unshift(@alts, $max_gt_conf_var->{'REF'});

              while($prev_gt =~ /([\/\|]*)(\d+)/g)
              {
                $new_gt_call .= $1.$alleles_hash{$alts[$2]};
              }

              $new_entry->{$sample}->{'GT'} = $new_gt_call;
            }
            else
            {
              $new_entry->{$sample}->{'GT'} = join("/", ('.') x $ploidy);
            }
          }

          # Set filter field
          vcf_add_filter_txt($new_entry, $multiallelic_tag);

          # Print new entry
          $vcf->print_entry($new_entry);
          $num_of_printed++;

          # Print merged entries
          for(my $var = $start; $var <= $end; $var++)
          {
            vcf_add_filter_txt($vars_at_same_pos[$var], $dupallele_tag);
            $vcf->print_entry($vars_at_same_pos[$var]);
          }

          my $vars_merged = $end-$start+1;
          $num_of_printed += $vars_merged;

          $biggest_merge = max($biggest_merge, $vars_merged);
          $num_of_vars_merged += $vars_merged;
        }

        # Update $i
        $i = $end;
      }
    }
  }
}

