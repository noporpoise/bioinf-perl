package IndexedString;

use strict;
use warnings;

use Carp;

sub new
{
  my $class = shift;
  my $kmer_size = shift;

  my $self = {'_ksize' => $kmer_size, '_strings' => [], '_kmers' => {}};
  bless($self, $class);

  $self->add_strings(@_);
  return $self;
}

sub add_strings
{
  my $self = shift;
  my $kmer_size = $self->{_ksize};

  for my $string (@_) {
    my $strindx = scalar(@{$self->{_strings}});
    push(@{$self->{_strings}}, $string);
    for(my $i = 0; $i < length($string)-$kmer_size+1; $i++) {
      my $kmer = substr($string, $i, $kmer_size);
      if(!defined($self->{_kmers}->{$kmer})) { $self->{_kmers}->{$kmer} = []; }
      push(@{$self->{_kmers}->{$kmer}}, [$strindx, $i]);
    }
  }
}

sub find_index
{
  my ($self, $str) = @_;
  my $kmer_size = $self->{_ksize};
  my $num_strings = scalar(@{$self->{_strings}});

  if(length($str) < $kmer_size)
  {
    my $idx;
    for(my $i = 0; $i < $num_strings; $i++) {
      if(($idx = index($self->{_strings}->[$i], $str)) > -1) {
        return ($i, $idx);
      }
    }
    return (-1,-1);
  }

  my $kmer0 = substr($str, 0, $kmer_size);
  my $kmers = $self->{_kmers}->{$kmer0};

  if(!defined($kmers)) { return (-1,-1); }

  for my $hit (@$kmers) {
    my ($idx,$pos) = @$hit;
    if(substr($self->{_strings}->[$idx], $pos, length($str)) eq $str) {
      return ($idx, $pos);
    }
  }

  return (-1,-1);
}

1;
