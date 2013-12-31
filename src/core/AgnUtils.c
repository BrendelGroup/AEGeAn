#include "AgnUtils.h"

GtArray* agn_array_copy(GtArray *source, size_t size)
{
  GtUword i;
  GtArray *new = gt_array_new(size);
  for(i = 0; i < gt_array_size(source); i++)
  {
    void *data = *(void **)gt_array_get(source, i);
    gt_array_add(new, data);
  }
  return new;
}

double agn_calc_splice_complexity(GtArray *transcripts)
{
  // FIXME
  return -1.0;
}

int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}