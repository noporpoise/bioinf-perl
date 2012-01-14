package VCFFile;

use strict;
use warnings;

use Carp;

# All methods are object methods except these:
use base 'Exporter';
our @EXPORT = qw(get_standard_vcf_columns);

my @header_tag_columns = qw(ALT FILTER FORMAT INFO);
my @header_tag_types = qw(Integer Float Character String Flag);
my @header_tag_hashkeys = qw(column ID Number Type Description);

sub new
{
  my $class = shift;
  my $handle = shift;
  my $next_line = <$handle>;
  
  if(!defined($next_line))
  {
    croak("VCF file is empty");
  }

  #
  # Load header
  #
  my %header_metainfo = ();
  my %header_tags = ();
  my @columns_arr = ();

  # Example header lines:
  ##fileDate=20090805
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001

  while(defined($next_line))
  {
    if($next_line =~ /##(.*)=(.*)/)
    {
      my ($key,$value) = ($1,$2);
      chomp($value);
      # header meta info line
      $header_metainfo{$key} = $value;
    }
    elsif($next_line =~ /##(\w+)=<(.*)>/)
    {
      my $tag;
      # Prints error and returns undef if not valid line
      if(defined($tag = _parse_header_tag($1,$2)))
      {
        $header_tags{$tag->{'ID'}} = $tag;
      }
    }
    elsif($next_line =~ /#[^#]/)
    {
      # column header line (e.g. '#CHROM..')
      my $header_line = substr($next_line, 1);
      chomp($header_line);

      @columns_arr = split(/\t+/, $header_line);

      if(@columns_arr == 0)
      {
        carp("VCF columns missing");
      }

      # Peak at first entry
      $next_line = <$handle>;
      last;
    }
    elsif($next_line =~ /^##/)
    {
      chomp($next_line);
      carp("VCF header line unrecognised '$next_line'");
    }
    elsif($next_line !~ /^\s*$/)
    {
      # Assume looking at first entry
      last;
    }

    $next_line = <$handle>;
  }

  if(@columns_arr == 0 && defined($next_line))
  {
    # No columns given, so assume some set of standard columns
    # Test this assumption - error if not true
    my @col_values = split(/\t/, $next_line);

    my @expected_cols = get_standard_vcf_columns();
    
    # Can be more columns that standard to include all samples,
    # but fewer needs to be reported (fatal)
    if(@col_values < @expected_cols)
    {
      croak("Invalid VCF - missing column headers and too few columns");
    }

    @columns_arr = @expected_cols;

    # Include samples in list of columns
    my $num_of_samples = scalar(@col_values) - scalar(@expected_cols);

    for(my $i = 1; $i <= $num_of_samples; $i++)
    {
      push(@columns_arr, 'sample'.$i);
    }
  }

  # Create a hash as another way to access headers
  my %columns_hash = ();

  for(my $i = 0; $i < @columns_arr; $i++)
  {
    $columns_hash{$columns_arr[$i]} = $i;
  }

  my $self = {
      _handle => $handle,
      _next_line => $next_line,
      _header_metainfo => \%header_metainfo,
      _header_tags => \%header_tags,
      _columns_hash => \%columns_hash,
      _columns_arr => \@columns_arr,
      _unread_entries => []
  };

  bless $self, $class;
  return $self;
}

sub _peak_line
{
  my ($self) = @_;
  return $self->{_next_line};
}

sub _read_line
{
  my ($self) = @_;
  my $temp_line = $self->{_next_line};
  my $handle = $self->{_handle};
  
  $self->{_next_line} = <$handle>;
  
  return $temp_line;
}

#
# Headers
#

sub _parse_header_tag
{
  # $column=<$tag_str>
  my ($tag_col,$tag_str) = @_;

  my %tag = ('column' => $tag_col);

  my @tag_parts = split(/\s*,\s*/, $tag_str);

  for my $tag_part (@tag_parts)
  {
    if($tag_part =~ /^ID\s*=\s*(.*)\s*$/i)
    {
      if(defined($tag{'ID'})) {
        carp("VCF header tag has multiple IDs");
      }
      $tag{'ID'} = $1;
    }
    elsif($tag_part =~ /^Number\s*=\s*(.*)\s*$/i)
    {
      if(defined($tag{'Number'})) {
        carp("VCF header tag has multiple Number values");
      }
      $tag{'Number'} = $1;
    }
    elsif($tag_part =~ /^Type\s*=\s*(.*)\s*$/i)
    {
      if(defined($tag{'Type'})) {
        carp("VCF header tag has multiple Type values");
      }
      $tag{'Type'} = $1;
    }
    elsif($tag_part =~ /^Description\s*=\s*\"(.*)\"\s*$/i)
    {
      if(defined($tag{'Description'})) {
        carp("VCF header tag has multiple Description values");
      }
      $tag{'Description'} = $1;
    }
  }

  # Check values - print error if needed
  return _check_valid_header_tag(\%tag) ? \%tag : undef;
}

# Returns 0 if invalid, 1 if valid
sub _check_valid_header_tag
{
  my ($tag) = @_;

  # Check all the things!

  # column
  if(!defined($tag->{'column'}))
  {
    die("VCF header tag 'column' not set");
  }
  elsif(!grep(/^$tag->{'column'}$/, @header_tag_columns))
  {
    carp("VCF header tag column not one of ".join(",", @header_tag_types)."\n");
    return 0;
  }

  # ID
  if(!defined($tag->{'ID'}))
  {
    carp("VCF header tag id missing");
    return 0;
  }
  elsif($tag->{'ID'} =~ /\s/)
  {
    carp("VCF header tag id contains whitespace characters '$tag->{'ID'}'\n");
    return 0;
  }

  if($tag->{'column'} !~ /^(?:ALT|FILTER)$/)
  {
    # Number
    if(!defined($tag->{'Number'}))
    {
      carp("VCF header tag 'Number' attribute is missing ('.' or an int plz)");
      return 0;
    }
    elsif($tag->{'Number'} !~ /^(?:\d+|\.)$/)
    {
      carp("VCF header tag number of arguments is not an +ve int " .
           "'$tag->{'Number'}'\n");
      return 0;
    }

    # Type
    if(!defined($tag->{'Type'}))
    {
      carp("VCF header tag 'Type' attribute is missing (e.g. @header_tag_types)");
      return 0;
    }
    elsif(!grep(/^$tag->{'Type'}$/, @header_tag_types))
    {
      carp("VCF header tag Type not one of ".join(",", @header_tag_types)."\n");
      return 0;
    }
  }
  elsif(defined($tag->{'Number'}) || defined($tag->{'Type'}))
  {
    carp("VCF header ALT/FILTER tags cannot have Number or Type attributes\n");
    return 0;
  }

  # Combinations
  if(defined($tag->{'Type'}) && $tag->{'Type'} eq "Flag" &&
     $tag->{'Number'} ne "0")
  {
    carp("VCF header type 'Flag' cannot have 'Number' other than 0");
    return 0;
  }

  return 1;
}

sub _cmp_header_tags
{
  my ($a, $b) = @_;

  my $cmp;
  
  for my $tag_field ('column','ID','Type','Description')
  {
    my $cmp = $a->{$tag_field} cmp $a->{$tag_field};
  
    if($cmp != 0)
    {
      return $cmp;
    }
  }

  return $a->{'Number'} <=> $b->{'Number'};
}

sub print_header
{
  my ($self, $out) = @_;

  # Open out handle to stdout, if not already defined
  if(!defined($out))
  {
    open($out, ">-");
  }

  # Print metainfo lines
  my $header_metainfo = $self->{_header_metainfo};

  for my $key (sort keys %$header_metainfo)
  {
    print $out "##$key=$header_metainfo->{$key}\n";
  }

  # Print tags
  my $header_tags = $self->{_header_tags};

  for my $key (sort _cmp_header_tags keys %$header_tags)
  {
    my $tag = $header_tags->{$key};

    if($tag->{'column'} =~ /^(?:ALT|FILTER)$/)
    {
      print $out "##" . $tag->{'column'} . "=<" .
                 "ID=" . $tag->{'ID'} . "," .
                 "Description=\"" . $tag->{'Description'} . "\">\n";
    }
    else
    {
      print $out "##" . $tag->{'column'} . "=<" .
                 "ID=" . $tag->{'ID'} . "," .
                 "Number=" . $tag->{'Number'} . "," .
                 "Type=" . $tag->{'Type'} . "," .
                 "Description=\"" . $tag->{'Description'} . "\">\n";
    }
  }

  # Print columns
  my @columns_arr = @{$self->{_columns_arr}};

  if(@columns_arr > 0)
  {
    print $out "#" . join("\t", @columns_arr) . "\n";
  }
}

sub get_header
{
  my ($self) = @_;

  my $header_str = "";

  # Print header to a string
  open(my $fh_str, '>', \$header_str)
    or die("Could not open string for writing");

  $self->print_header($fh_str);
  close($fh_str);

  return $header_str;
}

# Metainfo e.g.
##thing=value

sub add_header_metainfo
{
  my ($self, $key, $value) = @_;

  $self->{_header_metainfo}->{$key} = $value;
}

sub remove_header_metainfo
{
  my ($self, $key) = @_;

  delete($self->{_header_metainfo}->{$key});
}

sub get_header_metainfo
{
  my ($self) = @_;

  return %{$self->{_header_metainfo}};
}

# Tags e.g.
##INFO=<ID=AS,Number=1,Type=Float,Description="Woot">

sub add_header_tag
{
  my ($self, $tag_col, $tag_id, $tag_number, $tag_type, $tag_description)
    = @_;

  # INFO, FILTER, FORMAT.. column is in upper case
  $tag_col = uc($tag_col);

  # Integer, String.. lowercase with uppercase first letter
  $tag_type = lc($tag_type);
  substr($tag_type,0,1) = uc(substr($tag_type,0,1));

  my $tag = {'column' => $tag_col,
             'ID' => $tag_id,
             'Description' => $tag_description};

  if(defined($tag_number))
  {
    $tag->{'Number'} = $tag_number;
  }

  if(defined($tag_number))
  {
    $tag->{'Type'} = $tag_type;
  }

  if(_check_valid_header_tag($tag))
  {
    $self->{_header_tags}->{$tag_id} = $tag;
  }
}

sub remove_header_tag
{
  my ($self,$tag_id) = @_;

  delete($self->{_header_tags}->{$tag_id});
}

sub get_header_tags
{
  my ($self) = @_;
  return %{$self->{_header_tags}};
}

#
# Samples
#

sub get_list_of_sample_names
{
  my ($self) = @_;

  my @cols_array = $self->get_columns_array();

  my %usual_fields = ();
  my @standard_cols = get_standard_vcf_columns();

  for my $standard_col (@standard_cols) {
    $usual_fields{uc($standard_col)} = 1;
  }

  my @samples = grep {!defined($usual_fields{uc($_)})} @cols_array;
  return @samples;
}


#
# Columns
#
sub get_columns_array
{
  my ($self) = @_;
  return @{$self->{_columns_arr}};
}

sub get_columns_hash
{
  my ($self) = @_;
  return %{$self->{_columns_hash}};
}

sub set_columns_with_hash
{
  my ($self, $cols_hashref) = @_;
  
  my @cols_arr = sort {$cols_hashref->{$a} <=> $cols_hashref->{$b}}
                   keys %$cols_hashref;

  $self->{_columns_hash} = $cols_hashref;
  $self->{_columns_arr} = \@cols_arr;
}

sub set_columns_with_arr
{
  my ($self, $cols_arrref) = @_;

  my %cols_hash = ();

  for(my $i = 0; $i < @$cols_arrref; $i++)
  {
    $cols_hash{$cols_arrref->[$i]} = $i;
  }

  $self->{_columns_hash} = \%cols_hash;
  $self->{_columns_arr} = $cols_arrref;
}

# Static VCF method
sub get_standard_vcf_columns
{
  return qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
}

#
# Entries
#

sub unread_entry
{
  my ($self, $entry) = @_;

  push(@{$self->{_unread_entries}}, $entry);
}

sub read_entry
{
  my ($self) = @_;

  if(@{$self->{_unread_entries}} > 0)
  {
    return pop(@{$self->{_entry_buffered}});
  }

  my $vcf_line;

  # Read over empty line
  while(defined($vcf_line = $self->_read_line()) && $vcf_line =~ /^\s*$/) {}

  # no entries found
  if(!defined($vcf_line))
  {
    return undef;
  }

  chomp($vcf_line);

  my %entry = (); # store details in this hash
  my @entry_cols = split(/\t/, $vcf_line);

  my %vcf_columns = %{$self->{_columns_hash}};

  for my $col_name (keys %vcf_columns)
  {
    if($col_name ne "INFO")
    {
      $entry{$col_name} = @entry_cols[$vcf_columns{$col_name}];
    }
  }

  my $num_of_cols = scalar(@{$self->{_columns_arr}});

  if(@entry_cols < $num_of_cols)
  {
    croak("Not enough columns in VCF entry (ID: ".$entry{'ID'}."; " .
          "got " . @entry_cols . " columns, expected " . $num_of_cols . ")");
  }
  elsif(@entry_cols > $num_of_cols)
  {
    croak("Too many columns in VCF entry (ID: ".$entry{'ID'}."; " .
          "got " . @entry_cols . " columns, expected " . $num_of_cols . ")");
  }

  my %info_col = ();
  my @info_entries = split(";", $entry_cols[$vcf_columns{'INFO'}]);

  my %info_flags = ();

  for my $info_entry (@info_entries)
  {
    if($info_entry =~ /(.*)=(.*)/)
    {
      # key=value pair
      $info_col{$1} = $2;
    }
    else
    {
      # Flag
      $info_flags{$info_entry} = 1;
    }
  }

  $entry{'INFO'} = \%info_col;
  $entry{'INFO_flags'} = \%info_flags;

  # Auto-correct chromosome names
  if($entry{'CHROM'} !~ /^chr/)
  {
    if($entry{'CHROM'} =~ /^chr(.*)$/i) { # matches only with case-insensitive
      $entry{'CHROM'} = 'chr'.$1;
    }
    else {
      $entry{'CHROM'} = 'chr'.$entry{'CHROM'};
    }
  }

  # Correct SVLEN
  $entry{'INFO'}->{'SVLEN'} = length($entry{'ALT'}) - length($entry{'REF'});

  if(length($entry{'REF'}) != 1 || length($entry{'ALT'}) != 1)
  {
    # variant is not a SNP
    $entry{'true_REF'} = substr($entry{'REF'}, 1);
    $entry{'true_ALT'} = substr($entry{'ALT'}, 1);
    $entry{'true_POS'} = $entry{'POS'} + 1;
  }
  else
  {
    # SNP
    $entry{'true_REF'} = $entry{'REF'};
    $entry{'true_ALT'} = $entry{'ALT'};
    $entry{'true_POS'} = $entry{'POS'};
  }

  return \%entry;
}

sub print_entry
{
  my ($self, $entry, $out_handle) = @_;

  if(!defined($out_handle))
  {
    open($out_handle, ">-");
  }

  my @columns_arr = @{$self->{_columns_arr}};

  print $out_handle $entry->{$columns_arr[0]};

  for(my $i = 1; $i < @columns_arr; $i++)
  {
    if($columns_arr[$i] eq "INFO")
    {
      my $info_hashref = $entry->{'INFO'};
      my $flags_hashref = $entry->{'INFO_flags'};
      
      my @entries = map {$_ . "=" . $info_hashref->{$_}} keys %$info_hashref;
      push(@entries, keys %$flags_hashref);
      
      # Sort INFO entries
      @entries = sort {$a cmp $b} @entries;
      
      print $out_handle "\t" . join(";", @entries);
    }
    else
    {
      print $out_handle "\t" . $entry->{$columns_arr[$i]};
    }
  }

  print $out_handle "\n";
}

1;
