#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use NeedlemanWunsch;

#
# Read in fasta/fastq files, produce de Bruijn graph in dot format
#

#  maroon brown plum turquoise wheat
my @colours = qw(red green blue darkorange2 purple plum pink navy);
#my @colours = qw(red green blue);

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./resolve_debruijn.pl <kmer> <graph|call|both> [in.fa [in.fq ..]]

  Example:
    ./resolve_debruijn.pl 5 both graphs/test.fa > graphs/test.dot
    dot -Tpng graphs/test.dot > graphs/test.png\n";

  exit;
}

if(@ARGV < 2)
{
  print_usage();
}

my $kmer_size = shift;
my $action = lc(shift);
my @files = @ARGV;

if($kmer_size !~ /^\d+$/ || $kmer_size < 1)
{
  print_usage("Invalid kmer value '$kmer_size'");
}

my @actions = qw(graph call both);
if(!grep($action, @actions))
{
  print_usage("Action must be 'graph', 'call' or 'both'");
}

if(@files == 0)
{
  push(@files, "-");
}

# Create aligner
my $nw = new NeedlemanWunsch();

# colours uses the *indices* of the colours in the list above
# graph{$kmer} = {'prev_edges' => {}, 'next_edges' => {}, 'colours' => {}}
my %graph;

# Build graph
for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq,$qual);

  while((($title,$seq,$qual) = $fastn->read_next()) && defined($title))
  {
    load_read($seq);
  }

  close_fastn_file($fastn);
}

# Resolve structure
for my $file (@files)
{
  my $fastn = open_fastn_file($file);

  my ($title,$seq,$qual);

  while((($title,$seq,$qual) = $fastn->read_next()) && defined($title))
  {
    resolve_read($seq);
  }

  close_fastn_file($fastn);
}

if(grep($action, qw(graph both)))
{
  dump_graph();
}

if(grep($action, qw(call both)))
{
  invoke_bubble_caller();
}

$nw->destructor();




#
# Functions
#

sub load_read
{
  my ($read) = @_;

  if(length($read) < $kmer_size)
  {
    return;
  }

  for(my $i = 0; $i < length($read)-$kmer_size; $i++)
  {
    my $new_kmer = substr($read, $i, $kmer_size);

    if(!defined($graph{$new_kmer}))
    {
      $graph{$new_kmer} = {'prev' => {}, 'next' => {},
                           'colours' => {}, 'colends' => {},
                           'visited' => 0};
    }

    if($i > 0)
    {
      my $prev_base = substr($read, $i-1, 1);
      $graph{$new_kmer}->{'prev'}->{$prev_base} = 1;
    }

    if($i+$kmer_size < length($read))
    {
      my $next_base = substr($read, $i+$kmer_size, 1);
      $graph{$new_kmer}->{'next'}->{$next_base} = 1;
    }
  }
}

sub dump_graph
{
  print "digraph G {\n";

  # Print kmers
  for my $tmp_kmer (sort keys %graph)
  {
    # node properties
    my @properties = ();

    my @node_cols = sort keys %{$graph{$tmp_kmer}->{'colours'}};
    my @colends = sort keys %{$graph{$tmp_kmer}->{'colends'}};

    if(@colends > 0)
    {
      push(@properties, "shape=box");
    }

    if(@node_cols == @colours)
    {
      push(@properties, "style=filled", "bgcolor=gray");
    }

    print $tmp_kmer . (@properties == 0 ? '' : '['.join(',', @properties).']') . "\n";
  }

  # Print edges
  print "{ edge [dir=both, arrowtail=odot, arrowhead=dot]\n";

  for my $tmp_kmer (sort keys %graph)
  {
    my $next = substr($tmp_kmer, 1);
    for my $base (sort keys %{$graph{$tmp_kmer}->{'next'}})
    {
      my $nxt_kmer = $next . $base;
      # Print a coloured edge for each colour that is in both nodes
      my $num_printed = 0;
      for my $col (sort keys %{$graph{$tmp_kmer}->{'colours'}})
      {
        if(defined($graph{$nxt_kmer}->{'colours'}->{$col}))
        {
          print "$tmp_kmer -> $nxt_kmer [color=$colours[$col]];\n";
          $num_printed++;
        }
      }

      if($num_printed == 0)
      {
        print "$tmp_kmer -> $nxt_kmer [color=black];\n";
      }
    }
  }
  print "}\n";

  print "}\n";
}

sub resolve_read
{
  my ($read) = @_;
  
  # First: split read by recurring kmer
  # Loop through kmers
  my %kmers = ();

  for(my $i = 0; $i < length($read)-$kmer_size+1; $i++)
  {
    my $kmer = substr($read, $i, $kmer_size);

    if(defined(my $pos = $kmers{$kmer}))
    {
      # ...aaa....bbb.... split into:
      # ...aaa....bb      and
      #     aa....bbb....
      my $left = substr($read, 0, $i + $kmer_size - 1);
      my $right = substr($read, $pos+1);
      #print STDERR "$read -> $left $right [$kmer]\n";
      resolve_read($left);
      resolve_read($right);
      return;
    }
    else
    {
      $kmers{$kmer} = $i;
    }
  }

  # Made it through without any repeating kmers

  my @sole_junctions = ();
  my @pair_junctions = ();

  for(my $i = 1; $i < length($read) - $kmer_size; $i++)
  {
    my $kmer = substr($read, $i, $kmer_size);
    my $indegree = node_indegree($kmer);
    my $outdegree = node_outdegree($kmer);

    if($indegree > 1 && $outdegree > 1)
    {
      push(@sole_junctions, $i);
    }

    if($indegree > 1 || $outdegree > 1)
    {
      push(@pair_junctions, $i);
    }
  }

  for my $pos (@sole_junctions)
  {
    my $seq = substr($read, $pos-1, $kmer_size + 2);
    colour_nodes($seq);
  }

  if(@pair_junctions > 1)
  {
    # Resolvable pair_junctions
    # Between every pair of pair_junctions need to add colour
    for(my $i = 0; $i < @pair_junctions-1; $i++)
    {
      for(my $j = $i+1; $j < @pair_junctions; $j++)
      {
        my $start = $pair_junctions[$i] - 1;
        my $end = $pair_junctions[$j] + 1 + $kmer_size;
        my $seq = substr($read, $start, $end - $start);
        colour_nodes($seq);
      }
    }
  }
}

sub colour_nodes
{
  my ($seq) = @_;

  #print STDERR "$seq\n";

  my $num_kmers = length($seq) - $kmer_size + 1;

  # Get colours for nodes and all neighbours
  my @cols_used = ();
  my %colour_set = ();

  for(my $i = 0; $i < $num_kmers; $i++)
  {
    my %cols = ();

    my $kmer = substr($seq, $i, $kmer_size);
    my @neighbours = neighbouring_active_nodes($kmer);

    for my $neighbour (@neighbours)
    {
      for my $col (keys %{$graph{$neighbour}->{'colours'}})
      {
        $cols{$col} = 1;
      }
    }

    for my $col (keys %cols)
    {
      $colour_set{$col}++;
    }

    $cols_used[$i] = \%cols;
  }


  my $fkmer = substr($seq, 0, $kmer_size);
  my $lkmer = substr($seq, -$kmer_size);

  # Default to last colour

  # Check if this is redundant information
  for my $col (keys %colour_set)
  {
    my $count = $colour_set{$col};

    if($count == $num_kmers)
    {
      if(defined($graph{$fkmer}->{'colends'}->{$col}) &&
         defined($graph{$lkmer}->{'colends'}->{$col}))
      {
        # Dupe
        return;
      }
    }
  }

  # Pick a colour
  my $new_col = $#colours;

  if(keys %colour_set < @colours)
  {
    for(my $col = 0; $col < @colours; $col++)
    {
      if(!defined($colour_set{$col}))
      {
        # This is a colour we can use
        $new_col = $col;
        last;
      }
    }
  }

  # Nope -- we're first!

  for(my $i = 0; $i < $num_kmers; $i++)
  {
    my $kmer = substr($seq, $i, $kmer_size);

    if(defined($cols_used[$i]->{$new_col}))
    {
      # Set this node to not active
      set_all_cols($kmer);
    }
    else
    {
      $graph{$kmer}->{'colours'}->{$new_col} = 1;
    }
  }

  # Set first and last kmers
  $graph{$fkmer}->{'colends'}->{$new_col} = 1;
  $graph{$lkmer}->{'colends'}->{$new_col} = 1;
}

# Call het sites in a single sample
sub invoke_bubble_caller
{
  for my $kmer (keys %graph)
  {
    if(node_outdegree($kmer) == 2 && num_cols_set($kmer) < @colours)
    {
      bubble_caller($kmer);
    }
  }
}

sub bubble_caller
{
  my ($kmer) = @_;

  print STDERR "bubble_caller: $kmer\n";

  # by definition next_bases is 2 elements long
  my @next_bases = keys %{$graph{$kmer}->{'next'}};

  my @kmers = map {substr($kmer,1).$_} @next_bases;
  my @paths = map {$kmer.$_} @next_bases;
  my @extend = (1,1);
  my @col_hashes = ({},{});

  $graph{$kmers[0]}->{'visited'} = 1;
  $graph{$kmers[1]}->{'visited'} = 1;

  while($extend[0] || $extend[1])
  {
    for(my $i = 0; $i < 2; $i++)
    {
      if($extend[$i])
      {
        if(defined($kmers[$i] = next_kmer($kmers[$i], $col_hashes[$i])))
        {
          $paths[$i] .= substr($kmers[$i],-1);

          if($graph{$kmers[$i]}->{'visited'})
          {
            $extend[$i] = 0;
          }

          $graph{$kmers[$i]}->{'visited'} = 1;
        }
        else
        {
          $extend[$i] = 0;
        }
      }
    }
  }

  # Remove visited status
  reset_visited($paths[0]);
  reset_visited($paths[1]);

  # Check if last kmer is in other path
  my $lkmer0 = substr($paths[0], -$kmer_size);
  my $lkmer1 = substr($paths[1], -$kmer_size);

  my $pos;

  if(($pos = index($paths[1], $lkmer0)) != -1)
  {
    $paths[1] = substr($paths[1], 0, $pos + $kmer_size);
  }
  elsif(($pos = index($paths[0], $lkmer1)) != -1)
  {
    $paths[0] = substr($paths[0], 0, $pos + $kmer_size);
  }

  if($pos != -1)
  {
    #print STDERR "P1: $paths[0]\n";
    #print STDERR "P2: $paths[1]\n";
  
    my $out;
    open($out, ">&", \*STDERR) or die "Can't dup STDERR: $!";
    my $result = $nw->do_alignment(@paths);
    $nw->print_alignment($result, $out);
    close($out);
  }
}

# If no next kmer, returns undef
# If path splits, returns undef
# If next kmer has all colours set -- returns undef
sub next_kmer
{
  my ($kmer, $col_hash) = @_;

  my @cols = keys %$col_hash;

  my @next_bases = keys %{$graph{$kmer}->{'next'}};

  # Get set of potential next kmers that match colours
  my @next_kmers = ();

  for my $base (@next_bases)
  {
    my $tmp = substr($kmer,1) . $base;

    if(num_cols_set($tmp) < @colours)
    {
      my $col_match = 1;

      for my $col (@cols)
      {
        if(!defined($graph{$tmp}->{'colours'}->{$col}))
        {
          $col_match = 0;
          last;
        }
      }

      if($col_match)
      {
        push(@next_kmers, $tmp);
      }
    }
  }

  # Check if we have just one next_kmer
  # If so update colours and return it
  if(@next_kmers != 1)
  {
    return undef;
  }

  my $next_kmer = $next_kmers[0];

  for my $col (keys %{$graph{$next_kmer}->{'colends'}})
  {
    if(defined($col_hash->{$col}))
    {
      # Colour ending here
      delete($col_hash->{$col});
    }
    else
    {
      # Loop over prev nodes to see if we are starting on colour
      my $add = 1;

      # Colour shouldn't be in any prev nodes if we're to pick it up
      for my $prev_base (keys %{$graph{$next_kmer}->{'prev'}})
      {
        my $prev_kmer = $prev_base.substr($next_kmer, 0, -1);
        if(defined($graph{$prev_kmer}->{'colours'}->{$col}))
        {
          $add = 0;
          last;
        }
      }

      if($add)
      {
        # Colour staring here
        $col_hash->{$col} = 1;
        #print STDERR "$kmer Following colour $colours[$col]\n";
      }
    }
  }

  return $next_kmer;
}

sub reset_visited
{
  my ($path) = @_;
  for(my $i = 0; $i < length($path)-$kmer_size+1; $i++)
  {
    my $kmer = substr($path, $i, $kmer_size);
    $graph{$kmer}->{'visited'} = 0;
  }
}

sub node_indegree
{
  my ($kmer) = @_;
  return scalar keys %{$graph{$kmer}->{'prev'}};
}

sub node_outdegree
{
  my ($kmer) = @_;
  return scalar keys %{$graph{$kmer}->{'next'}};
}

sub num_cols_set
{
  my ($kmer) = @_;
  return scalar keys %{$graph{$kmer}->{'colours'}};
}

sub neighbouring_active_nodes
{
  my ($kmer) = @_;

  my @kmers = ();

  my @prev = map {$_ . substr($kmer, 0, -1)} keys %{$graph{$kmer}->{'prev'}};
  my @next = map {substr($kmer, 1) . $_} keys %{$graph{$kmer}->{'next'}};

  for my $kmer (@prev, @next)
  {
    if(num_cols_set($kmer) < @colours)
    {
      push(@kmers, $kmer);
    }
  }

  return @kmers;
}

sub set_all_cols
{
  my ($kmer) = @_;
  map {$graph{$kmer}->{'colours'}->{$_} = 1} 0..$#colours;
}

sub get_hash_of_kmers
{
  my ($read) = @_;
  my %kmers_hash = ();
  my @kmers = map {substr($read, $_, $kmer_size)} 0..(length($read)-$kmer_size-1);
  map {$kmers_hash{$_}++} @kmers;
  return (\%kmers_hash, \@kmers);
}
