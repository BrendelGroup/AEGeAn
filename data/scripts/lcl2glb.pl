#!/usr/local/env perl
use strict;

# lcl2glb.pl: convert glb2lcl.pl output from locus based coordinates to sequence
#             (chromosome or scaffold) based coordinates.
#
# Usage: perl lcl2glb.pl < locus-features-local.gff3 > locus-features-global.gff3

while(<STDIN>)
{
  if(m/^#/)
  {
    print unless(m/^##sequence-region/);
    next;
  }
  next if(m/\tlocus\t/);

  my @fields = split(/\t/);
  my ($ilocuspos) = $fields[8] =~ m/iLocus_pos=([^;]+)/;
  my ($ilseqid, $ilstart, $ilend) = $ilocuspos =~ m/(.+)_(\d+)-(\d+)/;
  my $offset = $ilstart - 1;
  $fields[0] = $ilseqid;
  $fields[3] += $offset;
  $fields[4] += $offset;
  $fields[8] =~ s/iLocus_pos=([^;]+);//;
  print join("\t", @fields);
}
