#include "AgnUtils.h"

double agn_calc_splice_complexity(GtArray *transcripts)
{
  // FIXME
  return -1.0;
}

int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}