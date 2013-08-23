#include "AgnCanonGeneStream.h"
#include "AgnGtExtensions.h"

struct AgnCanonGeneStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *cache;
};

#define canon_gene_stream_cast(GS)\
        gt_node_stream_cast(agn_canon_gene_stream_class(), GS)

static int canon_gene_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *error)
{
  AgnCanonGeneStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = canon_gene_stream_cast(ns);

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
      if(agn_gt_feature_node_is_gene_feature(current))
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

static void canon_gene_stream_free(GtNodeStream *ns)
{
  AgnCanonGeneStream *stream = canon_gene_stream_cast(ns);
  gt_node_stream_delete(stream->in_stream);
  gt_queue_delete(stream->cache);
}

const GtNodeStreamClass *agn_canon_gene_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnCanonGeneStream),
                                   canon_gene_stream_free,
                                   canon_gene_stream_next);
  }
  return nsc;
}

GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *ns;
  AgnCanonGeneStream *stream;
  gt_assert(in_stream);
  ns = gt_node_stream_create(agn_canon_gene_stream_class(), false);
  stream = canon_gene_stream_cast(ns);
  stream->in_stream = gt_node_stream_ref(in_stream);
  stream->cache = gt_queue_new();
  return ns;
}


