#include "AgnLocusFilterStream.h"
#include "AgnLocus.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnLocusFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *filters;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define locus_filter_stream_cast(GS)\
        gt_node_stream_cast(locus_filter_stream_class(), GS)

/**
 * @function Implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* locus_filter_stream_class(void);

/**
 * @function Class destructor.
 */
static void locus_filter_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they pass the provided filtering criteria.
 */
static int locus_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_locus_filter_stream_new(GtNodeStream *in_stream,
                                          GtArray *filters)
{
  GtNodeStream *ns;
  AgnLocusFilterStream *stream;
  gt_assert(in_stream);
  ns = gt_node_stream_create(locus_filter_stream_class(), false);
  stream = locus_filter_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);

  if(filters == NULL)
    stream->filters = gt_array_new( sizeof(AgnLocusFilter) );
  else
    stream->filters = gt_array_ref(filters);

  return ns;
}

static const GtNodeStreamClass *locus_filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnLocusFilterStream),
                                   locus_filter_stream_free,
                                   locus_filter_stream_next);
  }
  return nsc;
}

static void locus_filter_stream_free(GtNodeStream *ns)
{
  AgnLocusFilterStream *stream = locus_filter_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_array_delete(stream->filters);
}

static int locus_filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                    GtError *error)
{
  AgnLocusFilterStream *stream;
  GtFeatureNode *fn;
  GtUword i;
  gt_error_check(error);
  stream = locus_filter_stream_cast(ns);

  while(1)
  {
    bool keeplocus = true;
    int had_err = gt_node_stream_next(stream->in_stream, gn, error);
    if(had_err)
      return had_err;
    if(!*gn)
      return 0;

    fn = gt_feature_node_try_cast(*gn);
    if(!fn)
      return 0;

    gt_assert(gt_feature_node_has_type(fn, "locus"));
    for(i = 0; i < gt_array_size(stream->filters); i++)
    {
      AgnLocusFilter *filter = gt_array_get(stream->filters, i);
      if(agn_locus_filter_test(*gn, filter) == false)
      {
        keeplocus = false;
        break;
      }
    }

    if(!keeplocus)
    {
      gt_genome_node_delete(*gn);
      continue;
    }
  }

  return 0;
}

bool agn_locus_filter_stream_unit_test(GT_UNUSED AgnUnitTest *test)
{
  return false;
}
