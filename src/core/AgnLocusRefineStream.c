/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/queue_api.h"
#include "extended/sort_stream_api.h"
#include "AgnGeneStream.h"
#include "AgnLocus.h"
#include "AgnLocusStream.h"
#include "AgnLocusRefineStream.h"
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
  GtStr *source;
  GtUword count;
  GtQueue *locusqueue;
  GtGenomeNode *cache;
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
                                       GtArray *iloci, AgnLocus *orig);

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
 * @function FIXME
 */
static void
locus_refine_stream_mark_for_deletion(AgnLocusRefineStream *stream,
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

/**
 * @function FIXME
 */
static void locus_refine_stream_test_data(const char *filename, GtQueue *queue,
                                          GtUword delta);

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
  stream->source = gt_str_new_cstr("AEGeAn::AgnLocusStream");
  stream->count = 0;
  stream->locusqueue = gt_queue_new();
  stream->cache = NULL;
  return ns;
}

void agn_locus_refine_stream_set_idformat(AgnLocusRefineStream *stream,
                                          const char *format)
{
  agn_assert(stream && format);
  gt_str_delete(stream->idformat);
  stream->idformat = gt_str_new_cstr(format);
}

void agn_locus_refine_stream_set_source(AgnLocusRefineStream *stream,
                                        const char *source)
{
  agn_assert(stream && source);
  gt_str_delete(stream->source);
  stream->source = gt_str_new_cstr(source);
}

bool agn_locus_refine_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  locus_refine_stream_test_data("data/gff3/acep-syndrome.gff3", queue, 500);
  bool test1 = gt_queue_size(queue) == 3;
  bool test1a = false;
  if(test1)
  {
    GtGenomeNode *locus = gt_queue_get(queue);
    GtRange locusrange = gt_genome_node_get_range(locus);
    test1 = test1 && locusrange.start == 1915858 && locusrange.end == 1918866;
    const char *elen = gt_feature_node_get_attribute((GtFeatureNode *)locus,
                                                     "effective_length");
    if(elen)
      test1a = strcmp(elen, "11466") == 0;
    gt_genome_node_delete(locus);
    
    locus = gt_queue_get(queue);
    locusrange = gt_genome_node_get_range(locus);
    test1 = test1 && locusrange.start == 1916352 && locusrange.end == 1927323;
    gt_genome_node_delete(locus);
    
    locus = gt_queue_get(queue);
    locusrange = gt_genome_node_get_range(locus);
    test1 = test1 && locusrange.start == 1918157 && locusrange.end == 1921155;
    gt_genome_node_delete(locus);
  }
  agn_unit_test_result(test, "Atta cephalotes (syndrome): coords", test1);
  agn_unit_test_result(test, "Atta cephalotes (syndrome): elen", test1a);
  gt_queue_delete(queue);

  queue = gt_queue_new();
  locus_refine_stream_test_data("data/gff3/mrot-cst.gff3", queue, 500);
  bool test2 = gt_queue_size(queue) == 2;
  bool test2a = false;
  if(test2)
  {
    GtGenomeNode *locus = gt_queue_get(queue);
    GtRange locusrange = gt_genome_node_get_range(locus);
    test2 = test2 && locusrange.start == 9652 && locusrange.end == 19311;
    const char *elen = gt_feature_node_get_attribute((GtFeatureNode *)locus,
                                                     "effective_length");
    if(elen)
      test2a = strcmp(elen, "9660") == 0;
    gt_genome_node_delete(locus);
    
    locus = gt_queue_get(queue);
    locusrange = gt_genome_node_get_range(locus);
    test2 = test2 && locusrange.start == 11405 && locusrange.end == 18146;
    gt_genome_node_delete(locus);
  }
  agn_unit_test_result(test, "Megachile rotundata CST: coords", test2);
  agn_unit_test_result(test, "Megachile rotundata CST: elen", test2a);
  gt_queue_delete(queue);

  return agn_unit_test_success(test);
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

  GtUword numgenes = gt_array_size(bin);
  agn_assert(numgenes > 1);
  if(numgenes > 2)
    return false;

  gt_array_sort(bin, (GtCompare)agn_genome_node_compare);
  GtGenomeNode **gn1 = gt_array_get(bin, 0);
  GtGenomeNode **gn2 = gt_array_get(bin, 1);
  GtRange range1 = gt_genome_node_get_range(*gn1);
  GtRange range2 = gt_genome_node_get_range(*gn2);
  if(!gt_range_contains(&range1, &range2))
    return false;

  GtFeatureNode *fn1 = gt_feature_node_cast(*gn1);
  GtFeatureNode *fn2 = gt_feature_node_cast(*gn2);
  GtArray *exons = agn_typecheck_select(fn1, agn_typecheck_exon);
  if(gt_array_size(exons) <= 1)
  {
    gt_array_delete(exons);
    return false;
  }
  
  GtUword i;
  bool overlap = false;
  for(i = 0; i < gt_array_size(exons); i++)
  {
    GtGenomeNode **exon = gt_array_get(exons, i);
    GtRange exonrange = gt_genome_node_get_range(*exon);
    overlap = gt_range_overlap(&exonrange, &range2);
    if(overlap)
      break;
  }
  gt_array_delete(exons);
  if(overlap)
    return false;

  GtStr *seqid = gt_genome_node_get_seqid(*gn1);
  AgnLocus *locus = agn_locus_new(seqid);
  agn_locus_add_feature(locus, fn1);
  gt_genome_node_ref(*gn1);
  gt_array_add(iloci, locus);

  locus = agn_locus_new(seqid);
  agn_locus_add_feature(locus, fn2);
  gt_feature_node_add_attribute((GtFeatureNode *)locus, "intron_gene",
                                "true");
  gt_genome_node_ref(*gn2);
  gt_array_add(iloci, locus);

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
                                       GtArray *iloci, AgnLocus *orig)
{
  GtFeatureNode *origfn = gt_feature_node_cast(orig);
  GtRange origrange = gt_genome_node_get_range(orig);
  GtUword origro = 0;
  const char *rostr = gt_feature_node_get_attribute(origfn, "right_overlap");
  if(rostr)
    origro = atol(rostr);

  GtUword i;
  for(i = 0; i < gt_array_size(iloci); i++)
  {
    GtGenomeNode **gn = gt_array_get(iloci, i);
    GtFeatureNode *fn = gt_feature_node_cast(*gn);
    GtRange gnrange = gt_genome_node_get_range(*gn);
    if(origrange.start + stream->delta > gnrange.start)
      gnrange.start = origrange.start;
    else
      gnrange.start -= stream->delta;
    if(gnrange.end + stream->delta > origrange.end)
      gnrange.end = origrange.end;
    else
      gnrange.end += stream->delta;
    
    agn_locus_set_range(*gn, gnrange.start, gnrange.end);
    agn_assert(gt_range_contains(&origrange, &gnrange));
    if(i == 0)
    {
      char lenstr[32];
      sprintf(lenstr, "%lu", gt_range_length(&origrange) - origro);
      gt_feature_node_add_attribute(fn, "effective_length", lenstr);
    }
    gt_feature_node_set_source(fn, stream->source);
    gt_queue_add(stream->locusqueue, *gn);
  }
  return;
}

static void locus_refine_stream_free(GtNodeStream *ns)
{
  agn_assert(ns);
  AgnLocusRefineStream *stream = locus_refine_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_str_delete(stream->idformat);
  gt_str_delete(stream->source);
  gt_queue_delete(stream->locusqueue);
  if(stream->cache != NULL)
    gt_genome_node_delete(stream->cache);
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
  locus_refine_stream_extend(stream, iloci, gn);

  locus_refine_stream_mark_for_deletion(stream, gn);
  gt_array_delete(iloci);
  while(gt_array_size(bins) > 0)
  {
    GtArray **bin = gt_array_pop(bins);
    gt_array_delete(*bin);
  }
  gt_array_delete(bins);
  return 0;
}

static void
locus_refine_stream_mark_for_deletion(AgnLocusRefineStream *stream,
                                      GtGenomeNode *gn)
{
  if(stream->cache != NULL)
    gt_genome_node_delete(stream->cache);
  stream->cache = gn;
}

static void locus_refine_stream_mint(AgnLocusRefineStream *stream,
                                     AgnLocus *locus)
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
    agn_assert(gt_queue_size(stream->locusqueue) > 0);
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
    //gt_array_reverse(iloci);
  }
  return iloci;
}

static void locus_refine_stream_test_data(const char *filename, GtQueue *queue,
                                          GtUword delta)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();

  current_stream = gt_gff3_in_stream_new_unsorted(1, &filename);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = gt_sort_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtLogger *logger = gt_logger_new(true, "", stderr);
  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, delta);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_refine_stream_new(last_stream, delta, 1, true);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtError *error = gt_error_new();
  GtArray *loci = gt_array_new( sizeof(AgnLocus *) );
  current_stream = gt_array_out_stream_new(last_stream, loci, error);
  agn_assert(!gt_error_is_set(error));
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "error loading unit test data: %s\n", gt_error_get(error));
    exit(1);
  }

  gt_array_reverse(loci);
  while(gt_array_size(loci) > 0)
  {
    AgnLocus **locus = gt_array_pop(loci);
    gt_queue_add(queue, *locus);
  }

  while(gt_queue_size(streams) > 0)
  {
    GtNodeStream *ns = gt_queue_get(streams);
    gt_node_stream_delete(ns);
  }
  gt_queue_delete(streams);
  gt_logger_delete(logger);
  gt_error_delete(error);
  gt_array_delete(loci);
}
