/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/queue_api.h"
#include "AgnLocusRefineStream.h"
#include "AgnLocus.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

#define locus_refine_stream_cast(GS)\
        gt_node_stream_cast(locus_refine_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnLocusRefineStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtUword delta;
  GtUword minoverlap;
  bool by_cds;
  GtStr *idformat;
  GtUword count;
  GtQueue *locusqueue;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function FIXME
 */
static
GtArray *locus_refine_stream_bin_features(AgnLocusRefineStream *stream,
                                          GtFeatureNode *locus);

/**
 * @function FIXME
 */
static bool refine_locus_check_intron_genes(AgnLocusRefineStream *stream,
                                            GtArray *bin, GtArray *iloci);

/**
 * @function Implement the node stream interface.
 */
static const GtNodeStreamClass *locus_refine_stream_class(void);

/**
 * @function FIXME
 */
static void locus_refine_stream_extend(AgnLocusRefineStream *stream,
                                       GtArray *iloci);

/**
 * @function Destructor: release instance data.
 */
static void locus_refine_stream_free(GtNodeStream *ns);

/**
 * @function FIXME
 */
static int locus_refine_stream_handler(AgnLocusRefineStream *stream,
                                       GtGenomeNode *gn);

/**
 * @function Mint an ID for the given locus, tally counts of children and
 * grandchildren.
 */
static void locus_refine_stream_mint(AgnLocusRefineStream *stream, AgnLocus *locus);

/**
 * @function Feeds feature nodes of type ``locus`` to the output stream.
 */
static int locus_refine_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                    GtError *error);

/**
 * @function FIXME
 */
static GtArray *locus_refine_stream_resolve_bins(AgnLocusRefineStream *stream,
                                                 GtArray *bins);

//------------------------------------------------------------------------------
// Method definitions
//------------------------------------------------------------------------------

GtNodeStream *agn_locus_refine_stream_new(GtNodeStream *in_stream,
                                          GtUword delta, GtUword minoverlap,
                                          bool by_cds)
{
  GtNodeStream *ns = gt_node_stream_create(locus_refine_stream_class(), false);
  AgnLocusRefineStream *stream = locus_refine_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->delta = delta;
  stream->minoverlap = minoverlap;
  stream->by_cds = by_cds;
  stream->idformat = gt_str_new_cstr("locus%lu");
  stream->count = 0;
  stream->locusqueue = gt_queue_new();
  return ns;
}

void agn_locus_refine_stream_set_idformat(AgnLocusRefineStream *stream,
                                          const char *format)
{
  agn_assert(stream && format);
  gt_str_delete(stream->idformat);
  stream->idformat = gt_str_new_cstr(format);
  gt_queue_delete(stream->locusqueue);
}

bool agn_locus_refine_stream_unit_test(GT_UNUSED AgnUnitTest *test)
{
  return false;
}

static
GtArray *locus_refine_stream_bin_features(AgnLocusRefineStream *stream,
                                          GtFeatureNode *locus)
{
  GtArray *features = agn_feature_node_get_children(locus);
  GtUword numfeatures = gt_array_size(features);
  agn_assert(numfeatures >= 2);

  GtArray *bins = gt_array_new( sizeof(GtArray *) );
  GtArray *bin = gt_array_new( sizeof(GtGenomeNode *) );
  GtGenomeNode **g1 = gt_array_get(features, 0);
  gt_array_add(bin, *g1);
  gt_array_add(bins, bin);

  GtUword i, j, current_bin = 0;
  for(i = 1; i < numfeatures; i++)
  {
    GtGenomeNode **gn = gt_array_get(features, i);
    bin = *(GtArray **)gt_array_get(bins, current_bin);
    bool overlaps = false;
    for(j = 0; j < gt_array_size(bin); j++)
    {
      GtGenomeNode **test_gn = gt_array_get(bin, j);
      if(agn_overlap_ilocus(*gn, *test_gn, stream->minoverlap, stream->by_cds))
      {
        overlaps = true;
        gt_array_add(bin, *gn);
        break;
      }
    }
    if(!overlaps)
    {
      bin = gt_array_new( sizeof(GtFeatureNode *) );
      gt_array_add(bin, *gn);
      gt_array_add(bins, bin);
      current_bin++;
    }
  }
  gt_array_delete(features);
  return bins;
}

static bool refine_locus_check_intron_genes(AgnLocusRefineStream *stream,
                                            GtArray *bin, GtArray *iloci)
{
  agn_assert(bin);
  agn_assert(iloci);

  agn_assert(gt_array_size(bin) > 1);
  gt_array_sort(bin, (GtCompare)agn_genome_node_compare);
  GtGenomeNode **gn1 = gt_array_get(bin, 0);
  GtFeatureNode *fn1 = gt_feature_node_cast(*gn1);
  GtArray *exons = agn_typecheck_select(fn1, agn_typecheck_exon);
  if(gt_array_size(exons) <= 1)
  {
    gt_array_delete(exons);
    return false;
  }

  // First, see whether any intron from the first gene contains all other genes
  GtUword i, j;
  bool overlap = false;
  for(i = 0; i < gt_array_size(exons); i++)
  {
    GtGenomeNode **exon = gt_array_get(exons, i);
    GtRange exonrange = gt_genome_node_get_range(*exon);
    for(j = 1; j < gt_array_size(bin); j++)
    {
      GtGenomeNode **gn = gt_array_get(bin, j);
      GtRange feature_range = gt_genome_node_get_range(*gn);
      overlap = gt_range_overlap(&exonrange, &feature_range);
      if(overlap)
        break;
    }
    if(overlap == false)
      break;
  }
  gt_array_delete(exons);
  if(overlap)
    return false;

  // If so, next check for overlap between intron genes
  overlap = false;
  for(i = 1; i < gt_array_size(bin); i++)
  {
    GtGenomeNode **gn_i = gt_array_get(bin, i);
    GtRange range_i = gt_genome_node_get_range(*gn_i);
    for(j = i + 1; j < gt_array_size(bin); j++)
    {
      GtGenomeNode **gn_j = gt_array_get(bin, j);
      GtRange range_j = gt_genome_node_get_range(*gn_j);
      if(gt_range_overlap(&range_i, &range_j))
      {
        overlap = true;
        break;
      }
    }
    if(overlap)
      break;
  }
  if(overlap)
    return false;

  // If one gene has an intron containing all the other genes, and the genes
  // within the intron do not overlap, then we can create a distinct iLocus
  // for each gene.
  GtStr *seqid = gt_genome_node_get_seqid(*gn1);
  AgnLocus *locus = agn_locus_new(seqid);
  agn_locus_add_feature(locus, fn1);
  gt_array_add(iloci, locus);
  for(i = 1; i < gt_array_size(bin); i++)
  {
    GtGenomeNode **gn_i = gt_array_get(bin, i);
    GtFeatureNode *fn_i = gt_feature_node_cast(*gn_i);
    locus = agn_locus_new(seqid);
    agn_locus_add_feature(locus, fn_i);
    gt_feature_node_add_attribute((GtFeatureNode *)locus, "effective_length",
                                  "0");
    gt_array_add(iloci, locus);
  }
  return true;
}

static const GtNodeStreamClass *locus_refine_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnLocusRefineStream),
                                   locus_refine_stream_free,
                                   locus_refine_stream_next);
  }
  return nsc;
}

static void locus_refine_stream_extend(AgnLocusRefineStream *stream,
                                       GtArray *iloci)
{
  return;
}

static void locus_refine_stream_free(GtNodeStream *ns)
{
  agn_assert(ns);
  AgnLocusRefineStream *stream = locus_refine_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_str_delete(stream->idformat);
}

static int locus_refine_stream_handler(AgnLocusRefineStream *stream,
                                       GtGenomeNode *gn)
{
  GtFeatureNode *locus = gt_feature_node_cast(gn);
  if(gt_feature_node_number_of_children(locus) < 2)
  {
    gt_queue_add(stream->locusqueue, locus);
    return 0;
  }

  GtArray *bins = locus_refine_stream_bin_features(stream, locus);
  GtArray *iloci = locus_refine_stream_resolve_bins(stream, bins);
  locus_refine_stream_extend(stream, iloci);

  gt_genome_node_delete(gn);
  gt_array_delete(iloci);
  while(gt_array_size(bins) > 0)
  {
    GtArray **bin = gt_array_pop(bins);
    gt_array_delete(*bin);
  }
  gt_array_delete(bins);
  return 0;
}

static void locus_refine_stream_mint(AgnLocusRefineStream *stream, AgnLocus *locus)
{
  agn_assert(stream);
  stream->count++;

  char locusid[256];
  sprintf(locusid, gt_str_get(stream->idformat), stream->count);
  gt_feature_node_set_attribute((GtFeatureNode *)locus, "ID", locusid);

  GtArray *types = gt_array_new( sizeof(const char *) );
  GtHashmap *countsbytype = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                           gt_free_func);
  GtFeatureNode *feature;
  GtFeatureNode *locusfn = (GtFeatureNode *)locus;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(locusfn);
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    const char *childtype = gt_feature_node_get_type(feature);
    GtUword *num_of_type = gt_hashmap_get(countsbytype, childtype);
    if(num_of_type == NULL)
    {
      char *type = gt_cstr_dup(childtype);
      gt_array_add(types, type);
      num_of_type = gt_malloc( sizeof(GtUword) );
      (*num_of_type) = 0;
      gt_hashmap_add(countsbytype, type, num_of_type);
    }
    (*num_of_type)++;

    GtFeatureNodeIterator*
      subiter = gt_feature_node_iterator_new_direct(feature);
    GtFeatureNode *subfeature;
    for(subfeature = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature = gt_feature_node_iterator_next(subiter))
    {
      childtype = gt_feature_node_get_type(subfeature);
      num_of_type = gt_hashmap_get(countsbytype, childtype);
      if(num_of_type == NULL)
      {
        char *type = gt_cstr_dup(childtype);
        gt_array_add(types, type);
        num_of_type = gt_malloc( sizeof(GtUword) );
        (*num_of_type) = 0;
        gt_hashmap_add(countsbytype, type, num_of_type);
      }
      (*num_of_type)++;
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  GtUword i;
  for(i = 0; i < gt_array_size(types); i++)
  {
    const char **attrkey = gt_array_get(types, i);
    GtUword *attrvalue = gt_hashmap_get(countsbytype, *attrkey);
    char value[32];
    sprintf(value, "%lu", *attrvalue);
    gt_feature_node_set_attribute(locusfn, *attrkey, value);
  }

  gt_hashmap_delete(countsbytype);
  gt_array_delete(types);
}

static int locus_refine_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                    GtError *error)
{
  agn_assert(ns && gn && error);

  AgnLocusRefineStream *stream = locus_refine_stream_cast(ns);

  if(gt_queue_size(stream->locusqueue) > 0)
  {
    *gn = gt_queue_get(stream->locusqueue);
    locus_refine_stream_mint(stream, *gn);
    return 0;
  }

  int result = gt_node_stream_next(stream->in_stream, gn, error);
  if(result || !*gn)
    return result;

  if(gt_feature_node_try_cast(*gn))
  {
    locus_refine_stream_handler(stream, *gn);
    *gn = gt_queue_get(stream->locusqueue);
    locus_refine_stream_mint(stream, *gn);
  }

  return 0;
}

static GtArray *locus_refine_stream_resolve_bins(AgnLocusRefineStream *stream,
                                                 GtArray *bins)
{
  agn_assert(stream);
  agn_assert(bins);
  GtArray *iloci = gt_array_new( sizeof(AgnLocus *) );
  GtUword numbins = gt_array_size(bins);
  GtUword i;
  for(i = 0; i < numbins; i++)
  {
    GtArray **bin = gt_array_get(bins, i);
    GtUword numfeatures = gt_array_size(*bin);
    agn_assert(numfeatures > 0);
    if(numfeatures == 1)
    {
      GtGenomeNode **gn = gt_array_get(*bin, 0);
      GtFeatureNode *fn = gt_feature_node_cast(*gn);
      GtStr *seqid = gt_genome_node_get_seqid(*gn);
      AgnLocus *locus = agn_locus_new(seqid);
      agn_locus_add_feature(locus, fn);
      gt_genome_node_ref(*gn);  // Compensate for deletion of its former locus
      gt_array_add(iloci, locus);
    }
    else
    {
      if(refine_locus_check_intron_genes(stream, *bin, iloci) == false)
      {
        GtGenomeNode **gn = gt_array_get(*bin, 0);
        GtFeatureNode *fn = gt_feature_node_cast(*gn);
        GtStr *seqid = gt_genome_node_get_seqid(*gn);
        AgnLocus *locus = agn_locus_new(seqid);
        while(gt_array_size(*bin) > 0)
        {
          gn = gt_array_pop(*bin);
          fn = gt_feature_node_cast(*gn);
          agn_locus_add_feature(locus, fn);
          gt_genome_node_ref(*gn);  // Compensate for deletion of former parent
        }
        gt_array_add(iloci, locus);
      }
    }
  }
  if(gt_array_size(iloci) > 1)
  {
    gt_array_sort(iloci, (GtCompare)agn_genome_node_compare);
    gt_array_reverse(iloci);
  }
  return iloci;
}