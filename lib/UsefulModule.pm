package UsefulModule;

use strict;
use warnings;

use Carp;
use Scalar::Util qw(looks_like_number);

# Inherit from the "Exporter" module which handles exporting functions.
# Most procedural modules make use of this.

use base 'Exporter';

# When the module is invoked, export, by default, the function "hello" into 
# the namespace of the using code.

our @EXPORT = qw(num2str mem2str str2mem str2num round_int round_decimal
                 pretty_fraction binary_search_nearest trim
                 open_stdin);

=head1 NAME
 
num2str - Formats numbers so they are readable

=head1 SYNOPSIS

  use Num2Str;
  print num2str($num);
  print num2str($num, $separator);

=head1 DESCRIPTION

Format 

=head2 Functions

The following function is exported by default

=head3 hello

  print num2str(12443);
  print num2str(4523123.4521, " ");

Returns the number with a separator (',' by default) between every 3 digits.  
Example '12,443' and '4 523 123.4521'.  

=cut

sub num2str
{
  my $num = shift;
  my $sep = shift; # optional: separator between thousands [default: ',']
  my $num_decimals = shift; # optional up to this num decimals (default all)
  my $force_decimal = shift; # optional if != 0, force $num_decimals

  my @places = split(/\./, abs($num));
  my $intDigits = $places[0];

  # add commas or spaces to $intDigits

  my $initialPlaces = length($intDigits) % 3;
  
  if($initialPlaces == 0 && length($intDigits) >= 3)
  {
    $initialPlaces = 3;
  }
  
  my $strNum = substr($intDigits,0,$initialPlaces);
  
  for(my $i = $initialPlaces; $i <= length($intDigits)-3; $i += 3)
  {
    $strNum .= (defined($sep) ? $sep : ",") . substr($intDigits,$i,3);
  }
  
  if($num < 0)
  {
    $strNum = "-" . $strNum;
  }

  my $decStr = "";

  if(@places > 1)
  {
    if(!defined($num_decimals) || length($places[1]) <= $num_decimals)
    {
      $decStr = $places[1];
    }
    else
    {
      $decStr = substr($places[1], 0, $num_decimals);
    }
  }

  if(defined($force_decimal) && $force_decimal && length($decStr) < $num_decimals) {
    $decStr .= $decStr.("0"x($num_decimals-length($decStr)));
  }

  if(length($decStr) > 0) {
    $strNum .= ".".$decStr;
  }

  return $strNum;
}

sub mem2str
{
  my $num_bytes = shift;
  my $full_unit = shift;
  my $decimals = shift;

  $full_unit = defined($full_unit) && $full_unit;

  my @units = qw(B KB MB GB TB PB EB);
  my @full_units = qw(bytes kilobytes megabytes gigabytes terabytes
                      petabytes exabytes);

  my $unit;
  my $num_cpy = $num_bytes;
  for($unit = 0; $num_cpy >= 1024 && $unit < @units; $unit++)
  {
    $num_cpy /= 1024;
  }

  my $bytes_in_unit = 2**(10 * $unit);
  my $num_of_units = $num_bytes / $bytes_in_unit;

  return num2str($num_of_units, undef, $decimals) . " " .
         ($full_unit ? $full_units[$unit] : $units[$unit]);
}

# Parse 3M => 3*(1<<20)
sub str2mem
{
  my ($str) = @_;
  my $multiple = 1;
  my %units = ('K'=>1<<10,'M'=>1<<20,'G'=>1<<30,'T'=>1<<40,'P'=>1<<50,'E'=>1<<60);
  if($str =~ /(([EPTGMK])B?)$/i) {
    $multiple = $units{$2};
    $str = substr($str,0,-length($1));
  }
  if(!looks_like_number($str)) { die("Not numerical: '$str'"); }
  return $str * $multiple;
}

sub str2num
{
  my ($str) = @_;
  my $multiple = 1;
  my %units = ('K'=>10**3,'M'=>10**6,'G'=>10**9,'T'=>10**12,'P'=>10**15,'E'=>10**18);
  if($str =~ /([EPTGMK])$/i) {
    $multiple = $units{$1};
    $str = substr($str,0,-length($1));
  }
  if(!looks_like_number($str)) { die("Not numerical: '$str'"); }
  return $str * $multiple;
}

sub pretty_fraction
{
  my ($nominator, $denominator, $places) = @_;

  if(!defined($places))
  {
    $places = 2;
  }

  if($denominator == 0)
  {
    return num2str($nominator) . " / " . num2str($denominator);
  }

  my $percent = sprintf("%.".$places."f", 100 * $nominator / $denominator);

  return num2str($nominator) . " / " . num2str($denominator) . " " .
         "(" . $percent . "%)";
}

sub round_int
{
  my ($num, $round_to) = @_;

  if(!defined($round_to))
  {
    $round_to = 1;
  }

  return int($num / $round_to + 0.5) * $round_to;
}

sub round_decimal
{
  my ($num, $decimal_places) = @_;

  my $multiply = 10**$decimal_places;

  return int($num * $multiply + 0.5) / $multiply;
}

# Returns index of nearest value
sub binary_search_nearest
{
  # Note rounds up: 1.5 is 'nearer' to 2 than 1

  my ($arr, $search_value, $lower_bound, $upper_bound) = @_;

  if(@$arr == 0)
  {
    die("binary_search_nearest cannot an search empty array for value " .
        "'$search_value'");
  }
  elsif(!defined($lower_bound))
  {
    $lower_bound = 0;
    $upper_bound = scalar(@$arr) - 1;
  }

  #print "$lower_bound,$upper_bound,".@$arr."\n";

  if($lower_bound + 1 == $upper_bound)
  {
    # This bit makes it nearest rather than pure search
    my $dist_lower = $search_value - $arr->[$lower_bound];
    my $dist_upper = $arr->[$upper_bound] - $search_value;

    return ($dist_lower < $dist_upper ? $lower_bound : $upper_bound);
  }

  my $middle = int(($lower_bound + $upper_bound) / 2);

  if($search_value > $arr->[$middle])
  {
    return binary_search_nearest($arr, $search_value, $middle, $upper_bound);
  }
  elsif($search_value < $arr->[$middle])
  {
    return binary_search_nearest($arr, $search_value, $lower_bound, $middle);
  }
  else
  {
    return $middle;
  }
}

sub trim
{
  my ($str) = @_;

  $str =~ s/^\s+//;
  $str =~ s/\s+$//;

  return $str;
}

sub open_stdin
{
  my ($error_msg) = @_;

  my $stdin_handle;

  # -p checks STDIN is connected to a pipe
  if(!(-p STDIN) || !open($stdin_handle, "<&=STDIN"))
  {
    croak(defined($error_msg) ? $error_msg : "Cannot open STDIN");
  }

  return $stdin_handle;
}

1;
