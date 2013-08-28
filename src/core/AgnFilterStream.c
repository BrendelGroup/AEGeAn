#include "AgnFilterStream.h"
#include "AgnGtExtensions.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------

struct AgnFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *cache;
  GtHashmap *typestokeep;
  GtHashmap *typestofilter;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define filter_stream_cast(GS)\
        gt_node_stream_cast(filter_stream_class(), GS)

/**
 * Function that implements the GtNodeStream interface for this class.
 *
 * @returns    a node stream class object
 */
const GtNodeStreamClass* filter_stream_class(void);

/**
 * Destructor for the class.
 *
 * @param[in] ns    the node stream to be destroyed
 */
static void filter_stream_free(GtNodeStream *ns);

/**
 * Pulls nodes from the input stream and feeds them to the output stream if they
 * pass the provided filtering criteria.
 *
 * @param[in]  ns       the node stream
 * @param[out] gn       pointer to a genome node
 * @param[out] error    error object
 * @returns             0 in case of no error (*gn is set to next node or NULL
 *                      if stream is exhausted), -1 in case of error (error
 *                      object is set)
 */
static int filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream,
                                    GtHashmap *typestokeep,
                                    GtHashmap *typestofilter)
{
  GtNodeStream *ns;
  AgnFilterStream *stream;
  gt_assert(in_stream);
  gt_assert((typestokeep == NULL || typestofilter == NULL) &&
            (typestokeep != NULL || typestofilter != NULL));
  ns = gt_node_stream_create(filter_stream_class(), false);
  stream = filter_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->cache = gt_queue_new();
  stream->typestokeep = typestokeep;
  stream->typestofilter = typestofilter;
  return ns;
}

const GtNodeStreamClass *filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnFilterStream),
                                   filter_stream_free,
                                   filter_stream_next);
  }
  return nsc;
}

static void filter_stream_free(GtNodeStream *ns)
{
  AgnFilterStream *stream = filter_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_queue_delete(stream->cache);
  gt_hashmap_delete(stream->typestokeep);
  gt_hashmap_delete(stream->typestofilter);
}

static int filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error)
{
  AgnFilterStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = filter_stream_cast(ns);

  if(gt_queue_size(stream->cache) > 0)
  {
    *gn = gt_queue_get(stream->cache);
    return 0;
  }
  
  while(1)
  {
    had_err = gt_node_stream_next(stream->in_stream, gn, error);
    if(had_err)
      return had_err;
    if(!*gn)
      return 0;
    
    fn = gt_feature_node_try_cast(*gn);
    if(!fn)
      return 0;
    
    GtFeatureNode *current;
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      bool keepfeature;
      const char *type = gt_feature_node_get_type(current);
      if(stream->typestokeep == NULL)
      {
        keepfeature = true;
        if(gt_hashmap_get(stream->typestofilter, type) != NULL)
          keepfeature = false;
      }
      else
      {
        keepfeature = false;
        if(gt_hashmap_get(stream->typestokeep, type) != NULL)
          keepfeature = true;
      }

      if(keepfeature)
      {
        gt_genome_node_ref((GtGenomeNode *)current);
        gt_queue_add(stream->cache, current);
      }
    }
    gt_feature_node_iterator_delete(iter);
    gt_genome_node_delete((GtGenomeNode *)fn);
    if(gt_queue_size(stream->cache) > 0)
    {
      *gn = gt_queue_get(stream->cache);
      return 0;
    }
  }

  return 0;
}
