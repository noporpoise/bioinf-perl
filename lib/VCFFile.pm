package VCFFile;

use strict;
use warnings;

# All methods are object methods - nothing to export
use base 'Exporter';
our @EXPORT = qw(get_standard_vcf_columns vcf_add_to_header);

sub new
{
  my $class = shift;
  my $handle = shift;
  my $next_line = <$handle>;
  
  my $header = '';
  my %columns_hash = ();
  my @columns_arr = ();

  while(defined($next_line) && $next_line =~ /^##/)
  {
    $header .= $next_line;
    $next_line = <$handle>;
  }

  if($next_line !~ /#CHROM/)
  {
    if(length($header) > 0) {
      die("Invalid VCF - expected column headers");
    }
    else {
      # Assume some set of standard columns
      my @col_values = split(/\t/, $next_line);
      
      my @expected_cols = get_standard_vcf_columns();
      
      if(@col_values < @expected_cols) {
        die("Invalid VCF - missing column headers and too few columns");
      }

      @columns_arr = @expected_cols;
      
      for(my $i = 0; $i < @columns_arr; $i++) {
        $columns_hash{$columns_arr[$i]} = $i;
      }

      my $num_of_samples = scalar(@col_values) - scalar(@expected_cols);

      for(my $i = 1; $i <= $num_of_samples; $i++)
      {
        my $col_name = 'sample'.$i;
        push(@columns_arr, $col_name);
        $columns_hash{$col_name} = $col_name;
      }
    }
  }
  else
  {
    $header .= $next_line;
  
    my $header_line = substr($next_line, 1);
    chomp($header_line);
    
    @columns_arr = split(/\t/, $header_line);

    for(my $i = 0; $i < @columns_arr; $i++)
    {
      $columns_hash{$columns_arr[$i]} = $i;
    }

    # Peak at first entry
    $next_line = <$handle>;
  }
  
  my $self = {
      _handle => $handle,
      _next_line => $next_line,
      _header => $header,
      _columns_hash => \%columns_hash,
      _columns_arr => \@columns_arr
  };

  bless $self, $class;
  return $self;
}

sub peak_line
{
  my ($self) = @_;
  return $self->{_next_line};
}

sub read_line
{
  my ($self) = @_;
  my $temp_line = $self->{_next_line};
  my $handle = $self->{_handle};
  
  $self->{_next_line} = <$handle>;
  
  return $temp_line;
}

sub get_header
{
  my ($self) = @_;
  return $self->{_header};
}

sub get_columns
{
  my ($self) = @_;
  my %col_hash = %{$self->{_columns_hash}};
  my @col_arr = @{$self->{_columns_arr}};
  return (\%col_hash, \@col_arr);
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

sub read_entry
{
  my ($self) = @_;
  
  my $vcf_line = $self->read_line();
  
  if(!defined($vcf_line)) {
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
    die("Not enough columns in VCF entry " .
        "(got ".@entry_cols.", expected ".$num_of_cols.")");
  }
  elsif(@entry_cols > $num_of_cols)
  {
    die("Too many columns in VCF entry " .
        "(got ".@entry_cols.", expected ".$num_of_cols.")");
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
    if($entry{'CHROM'} =~ /^chr(.*)/i) { # matches only with case-insensitive
      $entry{'CHROM'} = 'chr'.$1;
    }
    else {
      $entry{'CHROM'} = 'chr'.$entry{'CHROM'};
    }
  }
  
  # OR use complex method in UsefulModule:
  #$entry{'CHROM'} = get_clean_chr_name($entry{'CHROM'});
  
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

  if(!defined($out_handle)) {
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

sub get_list_of_sample_names
{
  my ($self) = @_;

  my ($cols_hashref, $cols_arrref) = $self->get_columns();

  my %usual_fields = ();
  my @standard_cols = get_standard_vcf_columns();

  for my $standard_col (@standard_cols) {
    $usual_fields{uc($standard_col)} = 1;
  }

  my @samples = grep {!defined($usual_fields{uc($_)})} @$cols_arrref;
  return @samples;
}

sub get_standard_vcf_columns
{
  return qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
}

sub vcf_add_to_header
{
  my ($header,@new_lines) = @_;

  my @curr_lines = split(/[\n\r]+/, $header);

  if($curr_lines[$#curr_lines] !~ /#CHROM/i)
  {
    print STDERR "Warning: VCFFile.pm - " .
                 "header does not end with '#CHROM ..' line\n";
  }

  my $last_header_line = pop(@curr_lines);

  for my $new_line (@new_lines)
  {
    chomp($new_line);

    if($new_line !~ /^##/)
    {
      $new_line = '##'.$new_line;
      print STDERR "Warning: VCFFile.pm - " .
                   "adding hash characters to new header line\n";
    }

    push(@curr_lines, $new_line);
  }

  push(@curr_lines, $last_header_line);

  return join("\n", @curr_lines)."\n";
}

1;
