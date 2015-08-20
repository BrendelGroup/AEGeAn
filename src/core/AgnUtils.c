/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <string.h>
#include "core/hashmap_api.h"
#include "extended/feature_node_iterator_api.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"
#include "AgnVersion.h"

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
  // FIXME Function not yet (re)implemented
  return -1.0;
}

GtUword
agn_feature_index_copy_regions(GtFeatureIndex *dest, GtFeatureIndex *src,
                               bool use_orig, GtError *error)
{
  agn_assert(dest && src);
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
  agn_assert(dest && refrsrc && predsrc);
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
    gt_feature_index_has_seqid(predsrc, &predhas, seqid, error);

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
    else if(predhas)
      range_func(predsrc, &range, seqid, error);
    else
    {
      gt_str_delete(seqidstr);
      continue;
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

bool agn_feature_overlap_check(GtArray *feats)
{
  GtUword i,j;
  for(i = 0; i < gt_array_size(feats); i++)
  {
    GtFeatureNode *fn1 = *(GtFeatureNode **)gt_array_get(feats, i);
    GtRange range1 = gt_genome_node_get_range((GtGenomeNode *)fn1);
    for(j = i+1; j < gt_array_size(feats); j++)
    {
      GtFeatureNode *fn2 = *(GtFeatureNode **)gt_array_get(feats, j);
      GtRange range2 = gt_genome_node_get_range((GtGenomeNode *)fn2);
      if(gt_range_overlap(&range1, &range2))
        return true;
    }
  }
  return false;
}

GtRange agn_feature_node_get_cds_range(GtFeatureNode *fn)
{
  GtRange cds_range = {0,0};
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *child;
  for(child = gt_feature_node_iterator_next(iter);
      child != NULL;
      child = gt_feature_node_iterator_next(iter))
  {
    if(agn_typecheck_cds(child))
    {
      GtRange childrange = gt_genome_node_get_range((GtGenomeNode *)child);
      if(cds_range.end == 0)
        cds_range = childrange;
      else
        cds_range = gt_range_join(&cds_range, &childrange);
    }
  }
  gt_feature_node_iterator_delete(iter);
  return cds_range;
}

void agn_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)
{
  agn_assert(root && fn);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *child;
  for(child = gt_feature_node_iterator_next(iter);
      child != NULL;
      child = gt_feature_node_iterator_next(iter))
  {
    agn_feature_node_remove_tree(fn, child);
  }
  gt_feature_node_iterator_delete(iter);
  gt_feature_node_remove_leaf(root, fn);
}

int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

GtUword agn_mrna_3putr_length(GtFeatureNode *mrna)
{
  return agn_typecheck_feature_combined_length(mrna, agn_typecheck_utr3p);
}

GtUword agn_mrna_5putr_length(GtFeatureNode *mrna)
{
  return agn_typecheck_feature_combined_length(mrna, agn_typecheck_utr5p);
}

GtUword agn_mrna_cds_length(GtFeatureNode *mrna)
{
  return agn_typecheck_feature_combined_length(mrna, agn_typecheck_cds);
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

bool agn_overlap_ilocus(GtGenomeNode *f1, GtGenomeNode *f2,
                        GtUword minoverlap, bool by_cds)
{
  GtStr *seqid1 = gt_genome_node_get_seqid(f1);
  GtStr *seqid2 = gt_genome_node_get_seqid(f2);
  if(gt_str_cmp(seqid1, seqid2) != 0)
    return false;

  if(by_cds)
  {
    GtRange c1 = agn_feature_node_get_cds_range((GtFeatureNode *)f1);
    GtRange c2 = agn_feature_node_get_cds_range((GtFeatureNode *)f2);
    bool has_cds_1 = c1.end != 0;
    bool has_cds_2 = c2.end != 0;
    if(has_cds_1 != has_cds_2)
    {
      // One feature has a CDS, the other doesn't, so they should belong to
      // different iLoci even if they overlap.
      return false;
    }

    if(has_cds_1)
    {
      // Both have coding sequences, use those instead of the complete feature
      // coordinates.
      return gt_range_overlap_delta(&c1, &c2, minoverlap);
    }
  }

  // Either we are not in CDS mode, or the features don't have a CDS.
  GtRange r1 = gt_genome_node_get_range(f1);
  GtRange r2 = gt_genome_node_get_range(f2);
  return gt_range_overlap_delta(&r1, &r2, minoverlap);
}

void agn_print_version(const char *progname, FILE *outstream)
{
  fprintf(outstream, "[%s] AEGeAn Toolkit %s (%s %s)\n", progname,
          AGN_SEMANTIC_VERSION, AGN_VERSION_STABILITY, AGN_VERSION_HASH_SLUG);
}

int agn_sprintf_comma(GtUword n, char *buffer)
{
  if(n < 1000)
  {
    int spaces = sprintf(buffer, "%lu", n);
    buffer += spaces;
    return spaces;
  }
  int spaces = agn_sprintf_comma(n / 1000, buffer);
  buffer += spaces;
  sprintf(buffer, ",%03lu", n % 1000);
  return spaces + 4;
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
