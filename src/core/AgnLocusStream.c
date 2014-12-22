/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

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
  bool skip_empty;
  int endmode;
  GtFeatureIndex *seqranges;
  AgnLocus *prev_locus;
  GtQueue *locusqueue;
  GtGenomeNode *buffer;
  GtStr *source;
  GtStr *idformat;
  char *refrfile;
  char *predfile;
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
  stream->skip_empty = false;
  stream->endmode = 0;
  stream->seqranges = gt_feature_index_memory_new();
  stream->prev_locus = NULL;
  stream->locusqueue = gt_queue_new();
  stream->buffer = NULL;
  stream->source = gt_str_new_cstr("AEGeAn::AgnLocusStream");
  stream->idformat = gt_str_new_cstr("locus%lu");
  stream->refrfile = NULL;
  stream->predfile = NULL;
  return ns;
}

void agn_locus_stream_set_endmode(AgnLocusStream *stream, int endmode)
{
  agn_assert(stream);
  stream->endmode = endmode;
}

void agn_locus_stream_set_idformat(AgnLocusStream *stream, const char *format)
{
  agn_assert(stream && format);
  gt_str_delete(stream->idformat);
  stream->idformat = gt_str_new_cstr(format);
}

void agn_locus_stream_skip_empty_loci(AgnLocusStream *stream)
{
  agn_assert(stream);
  stream->skip_empty = true;
}

void agn_locus_stream_set_source(AgnLocusStream *stream, const char *source)
{
  agn_assert(stream && source);
  gt_str_delete(stream->source);
  stream->source = gt_str_new_cstr(source);
}

bool agn_locus_stream_unit_test(AgnUnitTest *test)
{
  return false;
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
  GtStr *prev_seqid;
  if(stream->prev_locus)
    prev_seqid = gt_genome_node_get_seqid(stream->prev_locus);

  // Handle initial loci
  if(stream->prev_locus == NULL || gt_str_cmp(seqid, prev_seqid) != 0)
  {
    if(locusrange.start >= seqrange.start + (2*stream->delta))
    {
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);
      if(stream->endmode >= 0 && !stream->skip_empty)
      {
        AgnLocus *emptylocus = agn_locus_new(seqid);
        agn_locus_set_range(emptylocus, seqrange.start,
                            locusrange.start - stream->delta - 1);
        gt_queue_add(stream->locusqueue, emptylocus);
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
    GtRange prev_range = gt_genome_node_get_range(stream->prev_locus);
    if(prev_range.end + stream->delta >= locusrange.start)
    {
      agn_locus_set_range(stream->prev_locus, prev_range.start,
                          locusrange.start - 1);
      agn_locus_set_range(locus, prev_range.end + 1, locusrange.end);
    }
    else if(prev_range.end + (2*stream->delta) >= locusrange.start)
    {
      agn_locus_set_range(stream->prev_locus, prev_range.start,
                          prev_range.end + stream->delta);
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);
    }
    else if(prev_range.end + (3*stream->delta) >= locusrange.start)
    {
      GtUword midpoint = (prev_range.end + locusrange.start) / 2;
      agn_locus_set_range(stream->prev_locus, prev_range.start, midpoint);
      agn_locus_set_range(locus, midpoint + 1, locusrange.end);
    }
    else
    {
      agn_locus_set_range(stream->prev_locus, prev_range.start,
                          prev_range.end + stream->delta);
      agn_locus_set_range(locus, locusrange.start - stream->delta,
                          locusrange.end);

      if(stream->endmode <= 0 && !stream->skip_empty)
      {
        AgnLocus *emptylocus = agn_locus_new(seqid);
        agn_locus_set_range(emptylocus,
                            prev_range.end + stream->delta + 1,
                            locusrange.start - stream->delta - 1);
        gt_queue_add(stream->locusqueue, emptylocus);
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
      if(stream->endmode >= 0)
      {
        AgnLocus *term_locus = agn_locus_new(seqid);
        agn_locus_set_range(term_locus,
                            locusrange.end + stream->delta + 1, seqrange.end);
        gt_queue_add(stream->locusqueue, term_locus);
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

    GtStr *newseqid = gt_genome_node_get_seqid(*gn);
    GtRange newrange = gt_genome_node_get_range(*gn);
    GtUword i;
    bool overlap = false;
    for(i = 0; i < gt_array_size(current_locus); i++)
    {
      GtGenomeNode **oldfeature = gt_array_get(current_locus, i);
      GtStr *oldseqid = gt_genome_node_get_seqid(*oldfeature);
      GtRange oldrange = gt_genome_node_get_range(*oldfeature);
      if(gt_str_cmp(newseqid, oldseqid) == 0 &&
         gt_range_overlap(&newrange, &oldrange))
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
  gt_str_delete(stream->idformat);
  gt_free(stream->refrfile);
  gt_free(stream->predfile);
}

static void locus_stream_mint(AgnLocusStream *stream, AgnLocus *locus)
{
  agn_assert(stream);
  stream->count++;

  char locusid[256];
  sprintf(locusid, gt_str_get(stream->idformat), stream->count);
  gt_feature_node_set_attribute((GtFeatureNode *)locus, "ID", locusid);
  gt_feature_node_set_source((GtFeatureNode *)locus, stream->source);

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

