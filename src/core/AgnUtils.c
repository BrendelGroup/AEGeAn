#include <string.h>
#include "core/hashmap_api.h"
#include "extended/feature_node_iterator_api.h"
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

GtUword
agn_feature_index_copy_regions(GtFeatureIndex *dest, GtFeatureIndex *src,
                               bool use_orig, GtError *error)
{
  gt_assert(dest && src);
  GtStrArray *seqids = gt_feature_index_get_seqids(src, error);
  GtUword i, rncount = 0;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtStr *seqidstr = gt_str_new_cstr(seqid);
    GtRange range;
    if(use_orig)
      gt_feature_index_get_orig_range_for_seqid(src, &range, seqid, error);
    else
      gt_feature_index_get_range_for_seqid(src, &range, seqid, error);

    GtGenomeNode *rn = gt_region_node_new(seqidstr, range.start, range.end);
    gt_feature_index_add_region_node(dest, (GtRegionNode *)rn, error);
    rncount++;

    gt_genome_node_delete(rn);
    gt_str_delete(seqidstr);
  }
  gt_str_array_delete(seqids);
  return rncount;
}

GtUword
agn_feature_index_copy_regions_pairwise(GtFeatureIndex *dest,
                                        GtFeatureIndex *refrsrc,
                                        GtFeatureIndex *predsrc,
                                        bool use_orig, GtError *error)
{
  gt_assert(dest && refrsrc && predsrc);
  GtStrArray *refr_seqids = gt_feature_index_get_seqids(refrsrc, error);
  GtStrArray *pred_seqids = gt_feature_index_get_seqids(predsrc, error);
  GtStrArray *seqids = agn_str_array_union(refr_seqids, pred_seqids);
  gt_str_array_delete(refr_seqids);
  gt_str_array_delete(pred_seqids);

  GtUword i, rncount = 0;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtStr *seqidstr = gt_str_new_cstr(seqid);
    GtRange range, refrrange, predrange;
    bool refrhas, predhas;
    gt_feature_index_has_seqid(refrsrc, &refrhas, seqid, error);
    gt_feature_index_has_seqid(refrsrc, &predhas, seqid, error);

    int (*range_func)(GtFeatureIndex *, GtRange *, const char *, GtError *);
    range_func = gt_feature_index_get_range_for_seqid;
    if(use_orig)
      range_func = gt_feature_index_get_orig_range_for_seqid;

    if(refrhas && predhas)
    {
      range_func(refrsrc, &refrrange, seqid, error);
      range_func(predsrc, &predrange, seqid, error);
      range = gt_range_join(&refrrange, &predrange);
    }
    else if(refrhas)
      range_func(refrsrc, &range, seqid, error);
    else
    {
      gt_assert(predhas);
      range_func(predsrc, &range, seqid, error);
    }

    GtGenomeNode *rn = gt_region_node_new(seqidstr, range.start, range.end);
    gt_feature_index_add_region_node(dest, (GtRegionNode *)rn, error);
    rncount++;

    gt_genome_node_delete(rn);
    gt_str_delete(seqidstr);
  }
  gt_str_array_delete(seqids);
  return rncount;
}

GtRange agn_multi_child_range(GtFeatureNode *top, GtFeatureNode *rep)
{
  GtRange range = {0, 0};
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(top);
  GtFeatureNode *child;
  for(child  = gt_feature_node_iterator_next(iter);
      child != NULL;
      child  = gt_feature_node_iterator_next(iter))
  {
    if(!gt_feature_node_is_multi(child))
      continue;

    if(gt_feature_node_get_multi_representative(child) == rep)
    {
      GtRange fnrange = gt_genome_node_get_range((GtGenomeNode *)child);
      if(range.start == 0 && range.end == 0)
        range = fnrange;
      else
        range = gt_range_join(&range, &fnrange);
    }
  }
  gt_feature_node_iterator_delete(iter);
  return range;
}

int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int agn_string_compare(const void *p1, const void *p2)
{
  const char *s1 = *(char **)p1;
  const char *s2 = *(char **)p2;
  return strcmp(s1, s2);
}

GtStrArray* agn_str_array_union(GtStrArray *a1, GtStrArray *a2)
{
  GtArray *strings = gt_array_new( sizeof(char *) );
  GtHashmap *added = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  GtUword i;
  for(i = 0; i < gt_str_array_size(a1); i++)
  {
    char *str = (char *)gt_str_array_get(a1, i);
    if(gt_hashmap_get(added, str) == NULL)
    {
      gt_hashmap_add(added, str, str);
      gt_array_add(strings, str);
    }
  }
  for(i = 0; i < gt_str_array_size(a2); i++)
  {
    char *str = (char *)gt_str_array_get(a2, i);
    if(gt_hashmap_get(added, str) == NULL)
    {
      gt_hashmap_add(added, str, str);
      gt_array_add(strings, str);
    }
  }
  gt_hashmap_delete(added);
  // The whole reason I'm going through this mess--GtStrArray class has no sort
  // function.
  gt_array_sort(strings, (GtCompare)agn_string_compare);

  GtStrArray *uniona = gt_str_array_new();
  for(i = 0; i < gt_array_size(strings); i++)
  {
    const char *str = *(const char **)gt_array_get(strings, i);
    gt_str_array_add_cstr(uniona, str);
  }
  gt_array_delete(strings);
  return uniona;
}
