#!/usr/local/env perl
use strict;

# glb2lcl.pl: convert LocusPocus output from sequence (chromosome or scaffold)
#             based coordinates to locus based coordinates.
#
# Usage: perl glb2lcl.pl < locus-features.gff3 > locus-features-local.gff3

my $ilocusid;
my $offset;
while(<STDIN>)
{
  if(m/^#/)
  {
    print unless(m/^##sequence-region/);
    next;
  }

  my @fields = split(/\t/);
  if($fields[2] eq "locus")
  {
    $ilocusid = sprintf("%s_%lu-%lu", $fields[0], $fields[3], $fields[4]);
    $offset = $fields[3] - 1;
    next;
  }
  $fields[0] = $ilocusid;
  $fields[3] -= $offset;
  $fields[4] -= $offset;
  $fields[8] =~ s/Parent=locus[^;]+;*//;
  print join("\t", @fields);
}
