package VCFFile;

use strict;
use warnings;

use List::Util qw(first min max);
use Carp;

# All methods are object methods except these:
use base 'Exporter';
our @EXPORT = qw(vcf_open vcf_get_standard_columns
                 vcf_sort_variants
                 vcf_add_filter_txt
                 vcf_is_snp vcf_get_clean_indel
                 vcf_get_ancestral_true_allele
                 vcf_get_flanks
                 vcf_get_alt_allele_genome_substr0
                 vcf_get_ancestral_genome_substr0
                 vcf_get_derived_genome_substr0
                 vcf_get_ref_alt_genome_lengths
                 vcf_get_slippage);

my @header_tag_columns = qw(ALT FILTER FORMAT INFO);
my @header_tag_types = qw(Integer Float Character String Flag);
my @header_tag_hashkeys = qw(column ID Number Type Description);

sub new
{
  my $class = shift;
  my $handle = shift;
  my $next_line = <$handle>;
  chomp($next_line);

  if(!defined($next_line))
  {
    # croak("VCF file is empty");
    $next_line = "#".join('\t', vcf_get_standard_columns())."\n";
  }

  #
  # Load header
  #
  # my @header_metainfo = ();
  # my %header_tags = ();
  my @header_extra_lines = (); # Unrecognised header lines
  my @columns_arr = ();

  # Example header lines:
  ##fileDate=20090805
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001

  while(defined($next_line))
  {
    if($next_line =~ /^##/)
    {
      my $tag = substr($next_line, 2);
      if($next_line =~ /^##(ALT|FILTER|INFO|FORMAT)=/i) {
        _check_valid_header_tag($tag);
      }
      push(@header_extra_lines, $tag);
    }
    elsif($next_line =~ /^#[^#]/)
    {
      # column header line (e.g. '#CHROM..')
      my $header_line = substr($next_line, 1);
      @columns_arr = split(/\t+/, $header_line);

      # Peek at first entry
      $next_line = <$handle>;

      if(@columns_arr == 0 || $columns_arr[0] ne "CHROM")
      {
        carp("VCF columns missing");
      }
      else { last; }
    }
    elsif($next_line !~ /^\s*$/)
    {
      # Assume looking at first entry
      last;
    }

    $next_line = <$handle>;
    chomp($next_line);
  }

  if(@columns_arr == 0 && defined($next_line))
  {
    # No columns given, so assume some set of standard columns
    # Test this assumption - error if not true
    my @col_values = split(/\t/, $next_line);

    my @expected_cols = vcf_get_standard_columns();

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

  # Get sample names
  my %usual_fields = ();
  my @standard_cols = vcf_get_standard_columns();

  for my $standard_col (@standard_cols) {
    $usual_fields{uc($standard_col)} = 1;
  }

  my @sample_names = grep {!defined($usual_fields{uc($_)})} @columns_arr;

  # Set _failed_vars_out to undef to skip non-PASS variants
  # Set _failed_vars_out to filehandle print non-PASS variants elsewhere
  # Do not set / delete() failed_vars_out to get all variants

  my $self = {
      _handle => $handle,
      _next_line => $next_line,
      _header_lines => \@header_extra_lines,
      _columns_hash => \%columns_hash,
      _columns_arr => \@columns_arr,
      _sample_names => \@sample_names,
      _unread_entries => []
  };

  bless $self, $class;
  return $self;
}

#
# Open VCF Handle
# If vcf_path is not defined, open stdin
#
sub vcf_open
{
  my ($vcf_path) = @_;
  my $vcffh;

  if(defined($vcf_path) && $vcf_path ne "-") {
    open($vcffh, $vcf_path) or croak("Cannot open VCF file '$vcf_path'\n");
  }
  elsif(-p STDIN) {
    # STDIN is connected to a pipe
    open($vcffh, "<&=STDIN") or croak("Cannot read pipe");
  }
  else
  {
    croak("Must specify or pipe in a VCF file");
  }

  my $vcf = new VCFFile($vcffh);
  return $vcf;
}

sub vcf_close
{
  my ($self) = @_;
  close($self->{_handle}) or carp("Cannot close VCF file handle");
}

sub _peek_line
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

sub set_filter_failed
{
  my ($self, $out_fh) = @_;

  $self->{_failed_vars_out} = $out_fh;
}

sub unset_filter_failed
{
  my ($self) = @_;
  delete($self->{_failed_vars_out});
}

#
# Headers
#

sub _parse_header_tag
{
  # $column=<$str>
  my ($tag_col, $str) = @_;
  $tag_col = uc($tag_col);

  if(!grep(/^$tag_col$/, @header_tag_columns) ||
     substr($str,0,1) ne '<' || substr($str,-1) ne '>')
  {
    carp("VCF header expected ##$tag_col=<...> but missing <>");
    return undef;
  }

  my %tag = ('column' => $tag_col);
  $str = substr($str, 1, -1);

  while($str =~ /\s*(\w+)\s*=\s*(\"(?:(?:\\\\)*\\\"|[^\"]*)*\"|\'(?:(?:\\\\)*\\\'|[^\']*)*\'|[^,]*?)(?:,|$)/gi)
  {
    my $key = lc($1);
    my $value = $2;

    if(substr($value,0,1) eq "'" || substr($value,0,1) eq '"')
    {
      $value = substr($value,1,-1);
      $value =~ s/\\\\/\\/g;
      $value =~ s/\\\"/\"/g;
      $value =~ s/\\\'/\'/g;
    }

    if($key eq "id")
    {
      if(defined($tag{'ID'})) {
        carp("VCF header tag has multiple IDs");
      }
      $tag{'ID'} = $value;
    }
    elsif($key eq "number")
    {
      if(defined($tag{'Number'})) {
        carp("VCF header tag has multiple Number values");
      }
      $tag{'Number'} = $value;
    }
    elsif($key eq "type")
    {
      if(defined($tag{'Type'})) {
        carp("VCF header tag has multiple Type values");
      }
      $tag{'Type'} = $value;
    }
    elsif($key eq "description")
    {
      if(defined($tag{'Description'})) {
        carp("VCF header tag has multiple Description values");
      }
      $tag{'Description'} = $value;
    }
    else
    {
      carp("VCF header tag has unknown key=value pair '$str'");
      return undef;
    }
  }

  return \%tag;
}

# Returns 0 if invalid, 1 if valid
sub _check_valid_header_tag
{
  my ($txt) = @_;

  my $tag;
  if($txt =~ /^(ALT|FILTER|INFO|FORMAT)=(.*)/i) {
    $tag = _parse_header_tag($1,$2);
  }

  if(!defined($tag)) { carp("Invalid header tag: $txt"); return 0; }

  # Check all the things!

  # ID
  if(!defined($tag->{'ID'}))
  {
    carp("VCF header tag id missing: '$txt'");
    return 0;
  }
  elsif($tag->{'ID'} =~ /\s/)
  {
    carp("VCF header tag id contains whitespace characters: '$txt'\n");
    return 0;
  }

  if($tag->{'column'} !~ /^(?:ALT|FILTER)$/)
  {
    # Number
    if(!defined($tag->{'Number'}))
    {
      carp("VCF header tag 'Number' attribute is missing ('.' or an int plz): '$txt'");
      return 0;
    }
    elsif($tag->{'Number'} !~ /^(?:\d+|\.)$/)
    {
      carp("VCF header tag Number= of arguments is not an integer or '.': '$txt'");
      return 0;
    }

    # Type
    if(!defined($tag->{'Type'}))
    {
      carp("VCF header tag 'Type' attribute is missing (e.g. @header_tag_types): '$txt'");
      return 0;
    }
    elsif($tag->{'column'} eq "INFO" &&
          !grep(/^$tag->{'Type'}$/, qw(Number Flag Character String)))
    {
      carp("VCF header tag Type not one of Number,Flag,Character,String: '$txt'");
      return 0;
    }
    elsif($tag->{'column'} eq "FORMAT" &&
          !grep(/^$tag->{'Type'}$/, qw(Integer Float Character String)))
    {
      carp("VCF header tag Type not one of Integer,Float,Character,String: '$txt'");
      return 0;
    }
  }
  elsif(defined($tag->{'Number'}) && $tag->{'Number'} != 0)
  {
    carp("VCF header ALT/FILTER tags cannot have Number attributes: '$txt'");
    return 0;
  }
  elsif(defined($tag->{'Type'}) && $tag->{'Type'} !~ /^FLAG$/i)
  {
    carp("VCF header ALT/FILTER tags cannot have Type attributes " .
         "('$tag->{'Type'}'): '$txt'");
    return 0;
  }

  if(!defined($tag->{'Description'}))
  {
    carp("VCF header tag missing Description (ID: $tag->{'ID'}): '$txt'");
    return 0;
  }

  # Combinations
  if(defined($tag->{'Type'}) && $tag->{'Type'} eq "Flag" &&
     $tag->{'Number'} ne "0")
  {
    carp("VCF header type 'Flag' cannot have 'Number' other than 0: '$txt'");
    return 0;
  }

  return 1;
}


sub print_header
{
  my ($self, $out) = @_;

  # Open out handle to stdout, if not already defined
  if(!defined($out)) { open($out, ">-"); }
 
  my $hdrs = $self->{_header_lines};

  my ($h1, $h2) = (-1,-1);

  if(defined($h1 = first {$hdrs->[$_] =~ /^fileformat=/i} 0..(@$hdrs-1))) {
    print "##".$hdrs->[$h1]."\n";
  }
  else { print "##fileformat=VCFv4.0\n"; $h1 = -1; }

  if(defined($h2 = first {$hdrs->[$_] =~ /^fileDate=/i} 0..(@$hdrs-1))) {
    print "##".$hdrs->[$h2]."\n";
  }
  else {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    print "##fileDate=".sprintf("%02d/%02d/%02d", $mday, $mon+1, $year % 100)."\n";
    $h2 = -1;
  }

  # Print in order of key e.g. contig=, ALT=, FORMAT=, INFO=, ...
  my %tags = ();
  for(my $i=0; $i<@$hdrs; $i++) {
    my ($tag) = ($hdrs->[$i] =~ /^[^=]*/);
    if(!defined($tags{$tag})) { $tags{$tag} = [$i]; }
    else { push(@{$tags{$tag}}, $i); }
  }

  for my $tag (sort {$a cmp $b} keys %tags) {
    for my $i (@{$tags{$tag}}) {
      if($i != $h1 && $i != $h2) { print "##".$hdrs->[$i]."\n"; }
    }
  }

  # Print header lines
  # for(my $i = 0; $i < @$hdrs; $i++) {
  #   if($i != $h1 && $i != $h2) { print "##".$hdrs->[$i]."\n"; }
  # }

  # Print columns
  my @columns_arr = @{$self->{_columns_arr}};
  if(@columns_arr > 0) {
    print $out "#" . join("\t", @columns_arr) . "\n";
  }
}

sub get_header
{
  my ($self) = @_;

  my $header_str = "";

  # Print header to a string
  my $fh_str;
  open($fh_str, '>', \$header_str)
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

  push(@{$self->{_header_lines}}, "$key=$value");
}

sub remove_header_metainfo
{
  my ($self, $key) = @_;

  my @newhdrs = grep {$_ !~ /^$key=/} @{$self->{_header_lines}};
  $self->{_header_lines} = \@newhdrs;
}

# Tags e.g.
##INFO=<ID=AS,Number=1,Type=Float,Description="Woot">

sub _prepare_header_tag
{
  my ($column, $tag_id, $tag_num, $tag_type, $description) = @_;

  # INFO, FILTER, FORMAT.. column is in upper case
  $column = uc($column);

  if(defined($tag_type)) {
    # Integer, String.. lowercase with uppercase first letter
    $tag_type = lc($tag_type);
    substr($tag_type,0,1) = uc(substr($tag_type,0,1));
  }

  $description =~ s/\\/\\\\/g;
  $description =~ s/\"/\\\"/g;

  if($column =~ /^(ALT|FILTER)$/i)
  {
    return "$column=<ID=$tag_id,Description=\"$description\">\n";
  }
  else
  {
    return "$column=<ID=$tag_id,Number=$tag_num,Type=$tag_type," .
           "Description=\"$description\">";
  }
}

sub add_header_tag
{
  my ($self, $tag_col, $tag_id, $tag_num, $tag_type, $tag_descr) = @_;

  my $tag = _prepare_header_tag($tag_col, $tag_id, $tag_num, $tag_type, $tag_descr);

  _check_valid_header_tag($tag);
  push(@{$self->{_header_lines}}, $tag);
}

sub remove_header_tag
{
  my ($self,$tagid) = @_;

  my $hdrs = $self->{_header_lines};
  my @newhdrs = grep {$_ !~ /^\w+=<(.*,)[\s]*ID=$tagid/} @$hdrs;
  $self->{_header_lines} = \@newhdrs;
}


#
# Samples
#

sub get_list_of_sample_names
{
  my ($self) = @_;
  return @{$self->{_sample_names}};
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
sub vcf_get_standard_columns
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

#
# my $vcf = new VCFFile("in.vcf");
# my $entry = $vcf->read_entry();
#
# $entry->{'CHROM'} is chromosome as in VCF
# $entry->{'POS'} is positions as in VCF
#
# $entry->{'REF'} is the reference allele as in the VCF
# $entry->{'ALT'} is the alternative allele as in the VCF
# $entry->{'true_REF'} is the bases of the ref allele that are actually affected
# $entry->{'true_ALT'} is the bases of the alt allele that are actually affected
#  e.g. REF='A';ALT='C';  true_REF='A';true_ALT='C' (SNP - no change)
#  e.g. REF='A';ALT='AC';  true_REF='';true_ALT='C' (indel)
#
# $entry->{'true_POS'} is the position (1-based) of the first affected base,
#   or the base AFTER a clean insertion.  In other words,
#   for SNPS $vcf_entry->{'POS'} == $vcf_entry->{'true_POS'}
#   for indels $vcf_entry->{'POS'} == $vcf_entry->{'true_POS'} - 1
#
#   $entry->{'true_POS'}+length($entry->{'true_REF'}) == base after the variant
#   $entry->{'true_POS'}-1 == base before the variant
#

sub read_entry
{
  my ($self) = @_;

  if(@{$self->{_unread_entries}} > 0)
  {
    return pop(@{$self->{_entry_buffered}});
  }

  my $entry;

  if(defined(my $fail_out = $self->{_failed_vars_out}))
  {
    # Print non-PASS variants to the given file handle
    while(defined($entry = $self->_read_entry_from_file()) &&
          defined($entry->{'FILTER'}) && $entry->{'FILTER'} ne "." &&
          uc($entry->{'FILTER'}) ne "PASS")
    {
      $self->print_entry($entry, $fail_out);
    }
  }
  elsif(exists($self->{_failed_vars_out}))
  {
    # Skip non-PASS variants
    while(defined($entry = $self->_read_entry_from_file()) &&
          defined($entry->{'FILTER'}) && $entry->{'FILTER'} ne "." &&
          uc($entry->{'FILTER'}) ne "PASS") {}
  }
  else
  {
    $entry = $self->_read_entry_from_file();
  }

  return $entry;
}

sub _read_entry_from_file
{
  my ($self) = @_;

  # store details in this hash
  my %entry = ();
  my @entry_cols;

  my $vcf_line;

  my %vcf_columns = %{$self->{_columns_hash}};

  # Read until we get the correct number of columns
  while(1)
  {
    # Read over empty line
    while(defined($vcf_line = $self->_read_line()) && $vcf_line =~ /^\s*$/) {}

    # no entries found
    if(!defined($vcf_line))
    {
      return undef;
    }

    chomp($vcf_line);

    # Reset entry
    %entry = ();
    @entry_cols = split(/\t/, $vcf_line);

    for my $col_name (grep {$_ ne "INFO"} keys %vcf_columns)
    {
      $entry{$col_name} = $entry_cols[$vcf_columns{$col_name}];
    }

    my $num_of_cols = scalar(@{$self->{_columns_arr}});

    if(@entry_cols < $num_of_cols)
    {
      warn("Not enough columns in VCF entry (ID: ".$entry{'ID'}."; " .
           "got " . @entry_cols . " columns, expected " . $num_of_cols . ")");
    }
    elsif(@entry_cols > $num_of_cols)
    {
      warn("Too many columns in VCF entry (ID: ".$entry{'ID'}."; " .
           "got " . @entry_cols . " columns, expected " . $num_of_cols . ")");
    }
    else
    {
      last;
    }
  }

  # Split up sample info
  my $samples_arr = $self->{_sample_names};
  my $format_str = $self->{_columns_hash}->{'FORMAT'};

  my @format_fields = split(":", $entry_cols[$vcf_columns{'FORMAT'}]);
  $entry{'FORMAT'} = \@format_fields;

  if(scalar(@$samples_arr) > 0 && defined($format_str))
  {
    for my $sample (@$samples_arr)
    {
      my @sample_fields = split(":", $entry{$sample});

      if(scalar(@sample_fields) != scalar(@format_fields))
      {
        croak("Sample does not match genotype format " .
              "[sample: '$sample'; var: " . $entry{'ID'} . "]");
      }

      my %entry_sample = ();
      @entry_sample{@format_fields} = @sample_fields;
      $entry{$sample} = \%entry_sample;
    }
  }

  # Get info data
  my %info_col = ();
  my %info_flags = ();

  my @info_entries = split(";", $entry_cols[$vcf_columns{'INFO'}]);

  if(!(@info_entries == 1 && $info_entries[0] eq "."))
  {
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
  }

  $entry{'INFO'} = \%info_col;
  $entry{'INFO_flags'} = \%info_flags;

  # Convert REF / ALT to upper case
  $entry{'REF'} = uc($entry{'REF'});
  $entry{'ALT'} = uc($entry{'ALT'});

  # Correct SVLEN
  $entry{'INFO'}->{'SVLEN'} = length($entry{'ALT'}) - length($entry{'REF'});

  # Correct AALEN
  if(defined($entry{'INFO'}->{'AA'}) && $entry{'INFO'}->{'AA'} =~ /^\d$/)
  {
    if($entry{'INFO'}->{'AA'} == 0)
    {
      $entry{'INFO'}->{'AALEN'} = length($entry{'ALT'}) - length($entry{'REF'});
    }
    elsif($entry{'INFO'}->{'AA'} == 1)
    {
      $entry{'INFO'}->{'AALEN'} = length($entry{'REF'}) - length($entry{'ALT'});
    }
    else
    {
      # Invalid AA value
      $entry{'INFO'}->{'AALEN'} = undef;
    }
  }

  # Strip padding bases off to get true alleles
  my $min_length = min(length($entry{'REF'}), length($entry{'ALT'}));
  my $padding_bases = 0;

  while($padding_bases < $min_length &&
        substr($entry{'REF'}, $padding_bases, 1) eq
          substr($entry{'ALT'}, $padding_bases, 1))
  {
    $padding_bases++;
  }

  if($padding_bases > 0)
  {
    # variant has a padding base
    $entry{'true_REF'} = substr($entry{'REF'}, $padding_bases);
    $entry{'true_ALT'} = substr($entry{'ALT'}, $padding_bases);
    $entry{'true_POS'} = $entry{'POS'} + $padding_bases;
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

# Print using columns read from VCF (set POS, REF, ALT - not true_X..)
sub print_entry
{
  my ($self, $entry, $out_handle) = @_;

  if(!defined($out_handle))
  {
    open($out_handle, ">-");
  }

  my @columns_arr = @{$self->{_columns_arr}};

  for(my $i = 0; $i < @columns_arr; $i++)
  {
    if($i > 0)
    {
      print $out_handle "\t";
    }

    if($columns_arr[$i] eq "INFO")
    {
      my $info_hashref = $entry->{'INFO'};
      my $flags_hashref = $entry->{'INFO_flags'};

      my @entries = map {$_ . "=" . $info_hashref->{$_}} keys %$info_hashref;
      push(@entries, keys %$flags_hashref);

      # Sort INFO entries
      @entries = sort {$a cmp $b} @entries;

      print $out_handle (@entries > 0 ? join(";", @entries) : '.');
    }
    elsif($columns_arr[$i] eq "FORMAT")
    {
      print $out_handle join(":", @{$entry->{'FORMAT'}});
    }
    elsif(ref($entry->{$columns_arr[$i]}) eq "HASH")
    {
      print $out_handle join(":", map {$entry->{$columns_arr[$i]}->{$_}}
                                      @{$entry->{'FORMAT'}});
    }
    else
    {
      print $out_handle $entry->{$columns_arr[$i]};
    }
  }

  print $out_handle "\n";
}

# sort variants (by chrom, pos, SVLEN,
#                   case insensitive ref-allele, case insensitive alt-allele)
sub cmp_variants
{
  my $order = $a->{'CHROM'} cmp $b->{'CHROM'};

  if(($order = $a->{'CHROM'} cmp $b->{'CHROM'}) != 0)
  {
    return $order;
  }

  if(($order = $a->{'POS'} <=> $b->{'POS'}) != 0)
  {
    return $order;
  }

  if(($order = $a->{'INFO'}->{'SVLEN'} <=> $b->{'INFO'}->{'SVLEN'}) != 0)
  {
    return $order;
  }

  if(($order = $a->{'REF'} cmp $b->{'REF'}) != 0)
  {
    return $order;
  }

  if(($order = $a->{'ALT'} cmp $b->{'ALT'}) != 0)
  {
    return $order;
  }

  return 0;
}

# sort variants (by chrom, pos, SVLEN,
#                   case insensitive ref-allele, case insensitive alt-allele)
sub vcf_sort_variants
{
  my ($variants) = @_;

  @$variants = sort cmp_variants @$variants;
}

# Add FILTER column txt
sub vcf_add_filter_txt
{
  my ($variant, $filter_txt) = @_;

  if($variant->{'FILTER'} eq "." || $variant->{'FILTER'} =~ /^PASS$/i)
  {
    $variant->{'FILTER'} = $filter_txt;
  }
  else
  {
    $variant->{'FILTER'} = $variant->{'FILTER'}.";".$filter_txt;
  }
}

sub get_ploidy
{
  my ($self, $var) = @_;

  for my $sample (@{$self->{_sample_names}})
  {
    if(defined($var->{$sample}->{'GT'}))
    {
      my @gts = split(/[\/\|]/, $var->{$sample}->{'GT'});
      return scalar(@gts);
    }
  }

  return undef;
}

# returns 0 or 1
sub vcf_is_snp
{
  my ($vcf_entry) = @_;

  my $ref_len = length($vcf_entry->{'true_REF'});
  my $alt_len = length($vcf_entry->{'true_ALT'});

  return ($ref_len == 1 && $alt_len == 1);
}

# returns undef or $indel
sub vcf_get_clean_indel
{
  my ($vcf_entry) = @_;

  my $ref = $vcf_entry->{'true_REF'};
  my $alt = $vcf_entry->{'true_ALT'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};

  if(min(length($ref), length($alt)) == 0 && $svlen != 0)
  {
    return $svlen > 0 ? $alt : $ref;
  }
  else
  {
    return undef;
  }
}

sub vcf_get_ancestral_true_allele
{
  my ($vcf_entry) = @_;

  my $ancestral = $vcf_entry->{'INFO'}->{'AA'};

  if(!defined($ancestral))
  {
    return undef;
  }

  return ($ancestral == 0 ? $vcf_entry->{'true_REF'} : $vcf_entry->{'true_ALT'});
}

sub vcf_get_flanks
{
  my ($vcf_entry, $genome, $flank_size) = @_;

  # Get position 0-based
  my $var_start = $vcf_entry->{'true_POS'} - 1;
  my $chr = $genome->guess_chrom_fasta_name($vcf_entry->{'CHROM'});

  if(!defined($chr))
  {
    return ("","");
  }

  my $ref_allele = $vcf_entry->{'true_REF'};

  my $left_flank_start = max(0, $var_start - $flank_size);
  my $left_flank_length = $var_start - $left_flank_start;

  my $right_flank_start = $var_start + length($ref_allele);
  my $right_flank_end = min($right_flank_start + $flank_size,
                            $genome->get_chr_length($chr));
  my $right_flank_length = $right_flank_end - $right_flank_start;

  my $left = $genome->get_chr_substr0($chr, $left_flank_start, $left_flank_length);
  my $right = $genome->get_chr_substr0($chr, $right_flank_start, $right_flank_length);

  return ($left, $right);
}

sub vcf_get_alt_allele_genome_substr0
{
  my ($vcf_entry, $chr, $genome, $start, $len) = @_;

  my $ref_start = $vcf_entry->{'true_POS'}-1;

  my $result = "";

  if($start < $ref_start)
  {
    my $up_to = min($start+$len, $ref_start);
    my $get_len = $up_to - $start;
    my $sub = $genome->get_chr_substr0($chr, $start, $get_len);

    $result .= $sub;
    $start += $get_len;
    $len -= $get_len;
  }

  if($len == 0)
  {
    return $result;
  }

  # Get derived allele sequence in place of reference
  my $alt_allele = $vcf_entry->{'true_ALT'};
  my $ref_len = length($vcf_entry->{'true_REF'});
  my $alt_len = length($alt_allele);

  if($start >= $ref_start && $start < $ref_start + $alt_len)
  {
    my $offset = $start-$ref_start;
    my $get_len = min($alt_len, $len) - $offset;
    my $sub = substr($alt_allele, $offset, $get_len);

    $result .= $sub;
    $start += $get_len;
    $len -= $get_len;
  }

  if($len == 0)
  {
    return $result;
  }

  # Get sequence after variant
  my $chr_len = $genome->get_chr_length($chr);
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};

  if($start >= $ref_start + $alt_len)
  {
    # Convert start to REF coordinates
    # svlen = length(alt) - length(ref)
    $start -= $svlen;

    my $up_to = min($start+$len, $chr_len);
    my $get_len = $up_to - $start;
    my $sub = $genome->get_chr_substr0($chr, $start, $get_len);

    $result .= $sub;
  }

  return $result;
}

sub vcf_get_ancestral_genome_substr0
{
  my ($vcf_entry, $chr, $genome, $start, $len) = @_;

  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  if($aa == 1)
  {
    return vcf_get_alt_allele_genome_substr0($vcf_entry, $chr, $genome, $start, $len);
  }
  else
  {
    return $genome->get_chr_substr0($chr, $start, $len);
  }
}

sub vcf_get_derived_genome_substr0
{
  my ($vcf_entry, $chr, $genome, $start, $len) = @_;

  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  if($aa == 0)
  {
    return vcf_get_alt_allele_genome_substr0($vcf_entry, $chr, $genome, $start, $len);
  }
  else
  {
    return $genome->get_chr_substr0($chr, $start, $len);
  }
}

sub vcf_get_ref_alt_genome_lengths
{
  my ($vcf_entry, $genome) = @_;
  my $chr_len = $genome->get_chr_length($vcf_entry->{'CHROM'});

  return  ($chr_len, $chr_len + $vcf_entry->{'INFO'}->{'SVLEN'});
}

# Uses left_flank, right_flank to get slippage (aka microhomology) for a variant
sub vcf_get_slippage
{
  my ($vcf_entry) = @_;

  my $long_allele = $vcf_entry->{'true_REF'};
  my $short_allele = $vcf_entry->{'true_ALT'};

  if(length($long_allele) < length($short_allele))
  {
    # swap
    ($long_allele, $short_allele) = ($short_allele, $long_allele);
  }

  my $left_flank_rev = reverse($vcf_entry->{'INFO'}->{'left_flank'});
  my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

  my $indel = $long_allele;
  my $indel_rev = reverse($indel);

  my $result = get_match($left_flank_rev, $indel_rev) ."".
               get_match($right_flank, $indel);

  return $result;
}

sub get_match
{
  my ($str1, $str2) = @_;

  $str1 = uc($str1);
  $str2 = uc($str2);

  my $len = min(length($str1), length($str2));

  my $i = 0;

  for($i = 0; $i < $len && substr($str1,$i,1) eq substr($str2,$i,1); $i++)
  {

  }

  return substr($str1, 0, $i);
}

1;
