/**

Copyright (c) 2010-2016, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <string.h>
#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "extended/sort_stream_api.h"
#include "AgnGeneStream.h"
#include "AgnInferParentStream.h"
#include "AgnLocusStream.h"
#include "AgnLocus.h"

#define locus_stream_cast(GS)\
        gt_node_stream_cast(locus_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnLocusStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtUword delta;
  GtUword count;
  bool skip_iiLoci;
  int endmode;
  GtFeatureIndex *seqranges;
  AgnLocus *prev_locus;
  GtQueue *locusqueue;
  GtGenomeNode *buffer;
  GtStr *source;
  GtStr *nameformat;
  char *refrfile;
  char *predfile;
  FILE *ilenfile;
};

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Add `feature` as a child of `locus`.
 */
static int locus_stream_add_feature(AgnLocusStream *stream, AgnLocus *locus,
                                    GtFeatureNode *feature, GtError *error);

/**
 * @function Implement the node stream interface.
 */
static const GtNodeStreamClass *locus_stream_class(void);

/**
 * @function Extend the locus coordinates.
 */
static void locus_stream_extend(AgnLocusStream *stream, AgnLocus *locus);

/**
 * @function Callback function: collect overlapping top-level features into
 * distinct interval loci.
 */
static int locus_stream_fn_handler(AgnLocusStream *stream, GtGenomeNode **gn,
                                   GtError *error);

/**
 * @function Destructor: release instance data.
 */
static void locus_stream_free(GtNodeStream *ns);

/**
 * @function Mint an ID for the given locus, tally counts of children and
 * grandchildren.
 */
static void locus_stream_mint(AgnLocusStream *stream, AgnLocus *locus);

/**
 * @function Feeds feature nodes of type ``locus`` to the output stream.
 */
static int locus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                             GtError *error);

/**
 * @function Callback function: store region nodes in a feature index to enable
 * computing end locus coordinates correctly.
 */
static int locus_stream_rn_handler(AgnLocusStream *stream, GtGenomeNode **gn,
                                   GtError *error);

/**
 * @function Load data from the following file(s) for unit testing.
 */
static void locus_stream_test_data(GtQueue *queue, int numfiles,
                                   const char **filenames, bool pairwise);

/**
 * @function Load data from the following file(s) for unit testing.
 */
static GtNodeStream *locus_stream_test_data2(GtFeatureIndex *iloci,
                                             GtNodeStream *ns, GtUword delta,
                                             bool skipends);

/**
 * @function Run unit tests for loci with delta > 0.
 */
static void locus_stream_unit_test_iloci(AgnUnitTest *test);

/**
 * @function Run unit tests for loci with delta=0.
 */
static void locus_stream_unit_test_loci(AgnUnitTest *test);

//------------------------------------------------------------------------------
// Method definitions
//------------------------------------------------------------------------------

void agn_locus_stream_label_pairwise(AgnLocusStream *stream,
                                     const char *refrfile, const char *predfile)
{
  agn_assert(stream);
  if(stream->refrfile != NULL)
    gt_free(stream->refrfile);
  stream->refrfile = gt_cstr_dup(refrfile);
  if(stream->predfile != NULL)
    gt_free(stream->predfile);
  stream->predfile = gt_cstr_dup(predfile);
}

GtNodeStream *agn_locus_stream_new(GtNodeStream *in_stream, GtUword delta)
{
  GtNodeStream *ns = gt_node_stream_create(locus_stream_class(), false);
  AgnLocusStream *stream = locus_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->delta = delta;
  stream->count = 0;
  stream->skip_iiLoci = false;
  stream->endmode = 0;
  stream->seqranges = gt_feature_index_memory_new();
  stream->prev_locus = NULL;
  stream->locusqueue = gt_queue_new();
  stream->buffer = NULL;
  stream->source = gt_str_new_cstr("AEGeAn::AgnLocusStream");
  stream->nameformat = NULL;
  stream->refrfile = NULL;
  stream->predfile = NULL;
  stream->ilenfile = NULL;
  return ns;
}

void agn_locus_stream_set_endmode(AgnLocusStream *stream, int endmode)
{
  agn_assert(stream);
  stream->endmode = endmode;
}

void agn_locus_stream_set_name_format(AgnLocusStream *stream, const char *fmt)
{
  agn_assert(stream && fmt);
  if(stream->nameformat)
    gt_str_delete(stream->nameformat);
  stream->nameformat = gt_str_new_cstr(fmt);
}

void agn_locus_stream_skip_iiLoci(AgnLocusStream *stream)
{
  agn_assert(stream);
  stream->skip_iiLoci = true;
}

void agn_locus_stream_set_source(AgnLocusStream *stream, const char *source)
{
  agn_assert(stream && source);
  gt_str_delete(stream->source);
  stream->source = gt_str_new_cstr(source);
}

bool agn_locus_stream_unit_test(AgnUnitTest *test)
{
  locus_stream_unit_test_loci(test);
  locus_stream_unit_test_iloci(test);
  return agn_unit_test_success(test);
}

void agn_locus_stream_track_ilens(AgnLocusStream *stream, FILE *ilenfile)
{
  stream->ilenfile = ilenfile;
}

static int locus_stream_add_feature(AgnLocusStream *stream, AgnLocus *locus,
                                    GtFeatureNode *feature, GtError *error)
{
  agn_assert(stream && locus && feature && error);
  if(stream->refrfile == NULL)
  {
    agn_locus_add_feature(locus, feature);
  }
  else
  {
    const char * filename = gt_genome_node_get_filename((GtGenomeNode*)feature);
    if(strcmp(filename, stream->refrfile) == 0)
      agn_locus_add_refr_feature(locus, feature);
    else if(strcmp(filename, stream->predfile) == 0)
      agn_locus_add_pred_feature(locus, feature);
    else
    {
      if(strcmp(filename, "generated") == 0)
      {
        gt_error_set(error, "cannot infer parent features while doing "
                     "comparative analysis; please preprocess annotations with "
                     "`canon-gff3` to create explicit `gene` features and then "
                     "try again");
      }
      else
      {
        gt_error_set(error, "filename '%s' does not match reference or "
                     "prediction", filename);
      }
      return -1;
    }
  }

  return 0;
}

static const GtNodeStreamClass *locus_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnLocusStream),
                                   locus_stream_free,
                                   locus_stream_next);
  }
  return nsc;
}

static void locus_stream_extend(AgnLocusStream *stream, AgnLocus *locus)
{
  agn_assert(stream && locus);
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange locusrange = gt_genome_node_get_range(locus);
  GtRange seqrange;
  gt_feature_index_get_range_for_seqid(stream->seqranges, &seqrange,
                                       gt_str_get(seqid), NULL);
  GtStr *prev_seqid = NULL;
  if(stream->prev_locus)
    prev_seqid = gt_genome_node_get_seqid(stream->prev_locus);

  // Handle initial loci
  if(stream->prev_locus == NULL || gt_str_cmp(seqid, prev_seqid) != 0)
  {
    if(locusrange.start >= seqrange.start + (2*stream->delta))
    {
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);
      if(stream->endmode >= 0 && !stream->skip_iiLoci)
      {
        AgnLocus *filocus = agn_locus_new(seqid);
        GtRange irange = { seqrange.start,
                           locusrange.start - stream->delta - 1 };
        agn_locus_set_range(filocus, irange.start, irange.end);
        gt_genome_node_add_user_data(filocus, "iiLocus_type",
                                     gt_cstr_dup("fiLocus"), gt_free_func);
        gt_queue_add(stream->locusqueue, filocus);
      }
    }
    else
    {
      agn_locus_set_range(locus, seqrange.start, locusrange.end);
    }
  }

  // Handle internal loci
  if(stream->prev_locus && gt_str_cmp(seqid, prev_seqid) == 0)
  {
    GtFeatureNode *prevfn = gt_feature_node_cast(stream->prev_locus);
    GtFeatureNode *locusfn = gt_feature_node_cast(locus);
    GtRange prev_range = gt_genome_node_get_range(stream->prev_locus);
    bool ovlp1 = prev_range.end + stream->delta >= locusrange.start;
    bool ovlp2 = prev_range.end >= locusrange.start;
    if(ovlp1 || ovlp2)
    {
      GtUword newend = prev_range.end + stream->delta;
      if(newend > seqrange.end)
        newend = seqrange.end;
      agn_locus_set_range(stream->prev_locus, prev_range.start, newend);

      GtUword newstart = locusrange.start - stream->delta;
      if(seqrange.start + stream->delta > locusrange.start)
        newstart = seqrange.start;
      agn_locus_set_range(locus, newstart, locusrange.end);

      GtUword overlap = newend - newstart + 1;
      char ovrlp[16];
      sprintf(ovrlp, "%lu", overlap);
      gt_feature_node_add_attribute(prevfn, "right_overlap", ovrlp);
      gt_feature_node_add_attribute(locusfn, "left_overlap", ovrlp);
      gt_feature_node_add_attribute(prevfn, "iiLocus_exception",
                                    "delta-overlap-gene");

      if(stream->ilenfile != NULL)
        fprintf(stream->ilenfile, "%s\t0\n", gt_str_get(seqid));
      gt_feature_node_add_attribute(prevfn, "riil", "0");
      gt_feature_node_add_attribute(locusfn, "liil", "0");
    }
    else if(prev_range.end + (2*stream->delta) >= locusrange.start)
    {
      GtUword newend = prev_range.end + stream->delta;
      if(newend > seqrange.end)
        newend = seqrange.end;
      agn_locus_set_range(stream->prev_locus, prev_range.start, newend);
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);

      GtUword overlap = prev_range.end-locusrange.start+(2*stream->delta)+1;
      char ovrlp[16];
      sprintf(ovrlp, "%lu", overlap);
      gt_feature_node_add_attribute(prevfn, "right_overlap", ovrlp);
      gt_feature_node_add_attribute(locusfn, "left_overlap", ovrlp);
      gt_feature_node_add_attribute(prevfn, "iiLocus_exception",
                                    "delta-overlap-delta");

      if(stream->ilenfile != NULL)
        fprintf(stream->ilenfile, "%s\t0\n", gt_str_get(seqid));
      gt_feature_node_add_attribute(prevfn, "riil", "0");
      gt_feature_node_add_attribute(locusfn, "liil", "0");
    }
    else if(prev_range.end + (3*stream->delta) >= locusrange.start)
    {
      GtUword midpoint = (prev_range.end + locusrange.start) / 2;
      agn_locus_set_range(stream->prev_locus, prev_range.start, midpoint);
      agn_locus_set_range(locus, midpoint + 1, locusrange.end);
      gt_feature_node_add_attribute(prevfn, "iiLocus_exception",
                                    "delta-re-extend");

      if(stream->ilenfile != NULL)
        fprintf(stream->ilenfile, "%s\t0\n", gt_str_get(seqid));
      gt_feature_node_add_attribute(prevfn, "riil", "0");
      gt_feature_node_add_attribute(locusfn, "liil", "0");
    }
    else
    {
      GtUword newend = prev_range.end + stream->delta;
      if(newend > seqrange.end)
        newend = seqrange.end;
      agn_locus_set_range(stream->prev_locus, prev_range.start, newend);
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);

      if(stream->endmode <= 0 && !stream->skip_iiLoci)
      {
        AgnLocus *iilocus = agn_locus_new(seqid);
        GtRange irange = { prev_range.end + stream->delta + 1,
                           locusrange.start - stream->delta - 1 };
        agn_locus_set_range(iilocus, irange.start, irange.end);

        if(stream->ilenfile != NULL)
        {
          fprintf(stream->ilenfile, "%s\t%lu\n", gt_str_get(seqid),
                  gt_range_length(&irange));
        }
        char iilocuslen[32];
        sprintf(iilocuslen, "%lu", gt_range_length(&irange));
        gt_feature_node_add_attribute(prevfn, "riil", iilocuslen);
        gt_feature_node_add_attribute(locusfn, "liil", iilocuslen);

        const char *orientstrs[] = { "FF", "FR", "RF", "RR" };
        int orient = agn_locus_inner_orientation(stream->prev_locus, locus);
        GtFeatureNode *iilocfn = gt_feature_node_cast(iilocus);
        gt_feature_node_set_attribute(iilocfn, "fg_orient", orientstrs[orient]);
        gt_queue_add(stream->locusqueue, iilocus);
      }
    }
  }
  gt_queue_add(stream->locusqueue, locus);

  // Handle terminal loci
  locusrange = gt_genome_node_get_range(locus);
  GtStr *buffer_seqid = NULL;
  if(stream->buffer)
    buffer_seqid = gt_genome_node_get_seqid(stream->buffer);
  if(stream->buffer == NULL || buffer_seqid == NULL ||
     gt_feature_node_try_cast(stream->buffer) == NULL ||
     gt_str_cmp(seqid, buffer_seqid) != 0)
  {
    if(seqrange.end > (2*stream->delta) &&
       locusrange.end <= seqrange.end - (2*stream->delta))
    {
      agn_locus_set_range(locus, locusrange.start,
                          locusrange.end + stream->delta);
      if(stream->endmode >= 0 && !stream->skip_iiLoci)
      {
        AgnLocus *filocus = agn_locus_new(seqid);
        agn_locus_set_range(filocus,
                            locusrange.end + stream->delta + 1, seqrange.end);
        gt_genome_node_add_user_data(filocus, "iLocus_type",
                                     gt_cstr_dup("fiLocus"), gt_free_func);
        gt_queue_add(stream->locusqueue, filocus);
      }
    }
    else
    {
      agn_locus_set_range(locus, locusrange.start, seqrange.end);
    }
  }
}

static int locus_stream_fn_handler(AgnLocusStream *stream, GtGenomeNode **gn,
                                   GtError *error)
{
  agn_assert(stream && gn && error);

  GtArray *current_locus = gt_array_new( sizeof(GtFeatureNode *) );
  if(stream->buffer != NULL)
  {
    gt_array_add(current_locus, stream->buffer);
    stream->buffer = NULL;
  }

  bool again = false;
  int result = 0;
  do
  {
    if(gt_array_size(current_locus) > 0)
      result = gt_node_stream_next(stream->in_stream, gn, error);
    if(!*gn || result)
    {
      if(!*gn)
        stream->buffer = NULL;
      break;
    }
    if(gt_feature_node_try_cast(*gn) == NULL)
    {
      stream->buffer = *gn;
      break;
    }

    GtUword i;
    bool overlap = false;
    for(i = 0; i < gt_array_size(current_locus); i++)
    {
      GtGenomeNode **oldfeature = gt_array_get(current_locus, i);
      if(agn_overlap_ilocus(*gn, *oldfeature, 1, false))
      {
        overlap = true;
        break;
      }
    }
    if(overlap || gt_array_size(current_locus) == 0)
    {
      gt_array_add(current_locus, *gn);
      again = true;
    }
    else
    {
      stream->buffer = *gn;
      again = false;
    }
  } while(again);
  if(result)
    return result;

  bool haderror = false;
  if(gt_array_size(current_locus) > 0)
  {
    GtGenomeNode **rep = gt_array_get(current_locus, 0);
    GtStr *seqid = gt_genome_node_get_seqid(*rep);
    AgnLocus *locus = agn_locus_new(seqid);
    while(gt_array_size(current_locus) > 0)
    {
      GtFeatureNode **fn = gt_array_pop(current_locus);
      if(locus_stream_add_feature(stream, locus, *fn, error) == -1)
      {
        haderror = true;
        break;
      }
    }

    if(stream->delta > 0)
      locus_stream_extend(stream, locus);

    stream->prev_locus = locus;
    if(gt_queue_size(stream->locusqueue) > 0)
      locus = gt_queue_get(stream->locusqueue);
    locus_stream_mint(stream, locus);
    *gn = locus;
  }
  gt_array_delete(current_locus);

  if(haderror)
    return -1;
  return 0;
}

static void locus_stream_free(GtNodeStream *ns)
{
  agn_assert(ns);
  AgnLocusStream *stream = locus_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_feature_index_delete(stream->seqranges);
  gt_queue_delete(stream->locusqueue);
  gt_str_delete(stream->source);
  if(stream->nameformat)
    gt_str_delete(stream->nameformat);
  gt_free(stream->refrfile);
  gt_free(stream->predfile);
}

static void locus_stream_mint(AgnLocusStream *stream, AgnLocus *locus)
{
  agn_assert(stream);
  stream->count++;

  gt_feature_node_set_source((GtFeatureNode *)locus, stream->source);
  if(stream->nameformat)
  {
    char locusname[256];
    sprintf(locusname, gt_str_get(stream->nameformat), stream->count);
    gt_feature_node_set_attribute((GtFeatureNode *)locus, "Name", locusname);
  }

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
    char key[128];
    char value[32];
    sprintf(key, "child_%s", *attrkey);
    sprintf(value, "%lu", *attrvalue);
    gt_feature_node_set_attribute(locusfn, key, value);
  }

  gt_hashmap_delete(countsbytype);
  gt_array_delete(types);
}

static int locus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                             GtError *error)
{
  agn_assert(ns && gn && error);

  AgnLocusStream *stream = locus_stream_cast(ns);
  if(gt_queue_size(stream->locusqueue) > 0)
  {
    *gn = gt_queue_get(stream->locusqueue);
    locus_stream_mint(stream, *gn);
    return 0;
  }

  if(stream->buffer != NULL)
  {
    if(gt_region_node_try_cast(stream->buffer))
    {
      *gn = stream->buffer;
      stream->buffer = NULL;
      stream->prev_locus = NULL;
      return locus_stream_rn_handler(stream, gn, error);
    }
    else if(gt_feature_node_try_cast(stream->buffer))
      return locus_stream_fn_handler(stream, gn, error);
    else
    {
      *gn = stream->buffer;
      stream->buffer = NULL;
      return 0;
    }
  }

  int result = gt_node_stream_next(stream->in_stream, gn, error);
  if(result || !*gn)
    return result;

  if(gt_feature_node_try_cast(*gn))
    return locus_stream_fn_handler(stream, gn, error);

  if(gt_region_node_try_cast(*gn))
    return locus_stream_rn_handler(stream, gn, error);

  return 0;
}

static int locus_stream_rn_handler(AgnLocusStream *stream, GtGenomeNode **gn,
                                   GtError *error)
{
  agn_assert(stream && gn && error);
  GtRegionNode *rn = gt_region_node_cast(*gn);
  return gt_feature_index_add_region_node(stream->seqranges, rn, error);
}

static void locus_stream_test_data(GtQueue *queue, int numfiles,
                                   const char **filenames, bool pairwise)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();

  current_stream = gt_gff3_in_stream_new_unsorted(numfiles, filenames);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = gt_sort_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtHashmap *type_parents = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                           gt_free_func);
  gt_hashmap_add(type_parents, gt_cstr_dup("mRNA"), gt_cstr_dup("gene"));
  current_stream = agn_infer_parent_stream_new(last_stream, type_parents);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtLogger *logger = gt_logger_new(true, "", stderr);
  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  if(pairwise)
  {
    agn_assert(numfiles == 2);
    agn_locus_stream_label_pairwise((AgnLocusStream *)current_stream,
                                    filenames[0], filenames[1]);
  }
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
  gt_hashmap_delete(type_parents);
}

static GtNodeStream *locus_stream_test_data2(GtFeatureIndex *iloci,
                                             GtNodeStream *ns, GtUword delta,
                                             bool skipends)
{
  agn_assert(iloci && ns);
  GtError *error = gt_error_new();

  GtNodeStream *lstream = agn_locus_stream_new(ns, delta);
  if(skipends)
    agn_locus_stream_set_endmode((AgnLocusStream *)lstream, -1);

  GtNodeStream *fstream = gt_feature_out_stream_new(lstream, iloci);
  gt_node_stream_pull(fstream, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[AgnIntervalLocusStream::interval_locus_stream_test_data] "
            "error processing node stream: %s\n", gt_error_get(error));
  }

  gt_node_stream_delete(lstream);
  gt_error_delete(error);
  return fstream;
}

static void locus_stream_unit_test_iloci(AgnUnitTest *test)
{
  GtFeatureIndex *iloci = gt_feature_index_memory_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  const char *infile = "data/gff3/ilocus.in.gff3";
  GtNodeStream *gff3 = gt_gff3_in_stream_new_unsorted(1, &infile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3);
  GtNodeStream *fstream = locus_stream_test_data2(iloci, gff3, 200, false);

  GtStrArray *seqids = gt_feature_index_get_seqids(iloci, error);
  if(gt_str_array_size(seqids) != 25)
  {
    agn_unit_test_result(test, "Delta=200", false);
    return;
  }

  const char *seqid = gt_str_array_get(seqids, 0);
  GtArray *seqloci = gt_feature_index_get_features_for_seqid(iloci,seqid,error);
  bool test1_1 = (gt_array_size(seqloci) == 1);
  if(test1_1)
  {
    AgnLocus *locus = *(AgnLocus **)gt_array_get(seqloci, 0);
    GtRange range = gt_genome_node_get_range(locus);
    test1_1 = range.start == 1 && range.end == 900;
  }
  agn_unit_test_result(test, "Delta=200", test1_1);
  gt_array_delete(seqloci);

  gt_feature_index_delete(iloci);
  gt_str_array_delete(seqids);
  gt_logger_delete(logger);
  gt_node_stream_delete(gff3);
  gt_node_stream_delete(fstream);
  gt_error_delete(error);
}

static void locus_stream_unit_test_loci(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();

  const char *filenames[] = { "data/gff3/grape-refr.gff3",
                              "data/gff3/grape-pred.gff3",
                              "" };
  locus_stream_test_data(queue, 2, filenames, true);
  GtUword starts[] = {    72, 10503, 22053, 26493, 30020, 37652, 42669, 48012,
                       49739, 55535, 67307, 77131, 83378, 88551 };
  GtUword ends[] =   {  5081, 11678, 23448, 29602, 33324, 38250, 45569, 48984,
                       54823, 61916, 69902, 81356, 86893, 92176 };
  GtUword numrefr[] = { 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1 };
  GtUword numpred[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1 };
  bool grapetest = gt_queue_size(queue) == 14;
  if(grapetest)
  {
    int i = 0;
    while(gt_queue_size(queue) > 0)
    {
      AgnLocus *locus = gt_queue_get(queue);
      GtRange range = gt_genome_node_get_range(locus);
      GtUword refrtrans = agn_locus_num_refr_mrnas(locus);
      GtUword predtrans = agn_locus_num_pred_mrnas(locus);
      grapetest = grapetest && starts[i] == range.start &&
                  ends[i] == range.end && numrefr[i] == refrtrans &&
                  numpred[i] == predtrans;
      i++;
      agn_locus_delete(locus);
    }
  }
  agn_unit_test_result(test, "grape test (pairwise)", grapetest);


  filenames[0] = "data/gff3/pd0159-refr.gff3";
  filenames[1] = "data/gff3/pd0159-pred.gff3";
  locus_stream_test_data(queue, 2, filenames, true);
  GtUword pdstarts[] = { 15005, 25101, 27822,  33635,  40258,  42504, 50007,
                         56261, 60860, 73343,  93338, 107687, 107919 };
  GtUword pdends[]   = { 24351, 25152, 29494,  38145,  42162,  45986, 51764,
                         59660, 69505, 90631, 107441, 107862, 111581 };
  GtUword pdnumrefr[] = { 1, 0, 1, 0, 1, 1, 1, 1, 3, 1, 1, 0, 1 };
  GtUword pdnumpred[] = { 2, 1, 1, 1, 0, 1, 1, 1, 3, 3, 2, 1, 1 };
  bool pdtest = gt_queue_size(queue) == 13;
  if(pdtest)
  {
    int i = 0;
    while(gt_queue_size(queue) > 0)
    {
      AgnLocus *locus = gt_queue_get(queue);
      GtRange range = gt_genome_node_get_range(locus);
      GtUword refrtrans = agn_locus_num_refr_mrnas(locus);
      GtUword predtrans = agn_locus_num_pred_mrnas(locus);
      pdtest = pdtest && pdstarts[i] == range.start && pdends[i] == range.end &&
               pdnumrefr[i] == refrtrans && pdnumpred[i] == predtrans;
      i++;
      agn_locus_delete(locus);
    }
  }
  agn_unit_test_result(test, "Pdom test (pairwise)", pdtest);


  filenames[0] = "data/gff3/amel-aug-nvit-param.gff3",
  filenames[1] = "data/gff3/amel-aug-dmel-param.gff3",
  filenames[2] = "data/gff3/amel-aug-athal-param.gff3";
  locus_stream_test_data(queue, 3, filenames, false);
  GtUword augstarts[] = {     1, 36466, 44388, 72127, 76794 };
  GtUword augends[]   = { 33764, 41748, 70877, 76431, 97981 };
  GtUword augntrans[] = { 6, 3, 4, 2, 6 };
  bool augtest = gt_queue_size(queue) == 5;
  if(augtest)
  {
    int i = 0;
    while(gt_queue_size(queue) > 0)
    {
      AgnLocus *locus = gt_queue_get(queue);
      GtRange range = gt_genome_node_get_range(locus);
      GtUword ntrans = agn_locus_num_mrnas(locus);
      augtest = augtest &&
                augstarts[i] == range.start && augends[i] == range.end &&
                augntrans[i] == ntrans;
      i++;
      agn_locus_delete(locus);
    }
  }
  agn_unit_test_result(test, "Amel test (Augustus)", augtest);

  gt_queue_delete(queue);
}
