#include "AgnCanonGeneStream.h"
#include "AgnFilterStream.h"
#include "AgnGtExtensions.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnCanonGeneStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *streams;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

#define canon_gene_stream_cast(GS)\
        gt_node_stream_cast(canon_gene_stream_class(), GS)

/**
 * Function that implements the GtNodeStream interface for this class.
 *
 * @returns    a node stream class object
 */
static const GtNodeStreamClass* canon_gene_stream_class(void);

/**
 * Destructor for the class.
 *
 * @param[in] ns    the node stream to be destroyed
 */
static void canon_gene_stream_free(GtNodeStream *ns);

/**
 * Pulls nodes from the input stream and feeds them to the output stream if they
 * are valid canonical protein-coding genes.
 *
 * @param[in]  ns       the node stream
 * @param[out] gn       pointer to a genome node
 * @param[out] error    error object
 * @returns             0 in case of no error (*gn is set to next node or NULL
 *                      if stream is exhausted), -1 in case of error (error
 *                      object is set)
 */
static int canon_gene_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream,
                                        AgnLogger *logger)
{
  GtNodeStream *ns;
  AgnCanonGeneStream *stream;
  gt_assert(in_stream);

  ns = gt_node_stream_create(canon_gene_stream_class(), false);
  stream = canon_gene_stream_cast(ns);
  stream->streams = gt_queue_new();

  GtNodeVisitor *icnv = agn_infer_cds_visitor_new(logger);
  GtNodeVisitor *ienv = agn_infer_exons_visitor_new(logger);
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "gene", "gene");

  gt_queue_add(stream->streams, gt_node_stream_ref(in_stream));
  GtNodeStream *icnv_stream = gt_visitor_stream_new(in_stream, icnv);
  gt_queue_add(stream->streams, icnv_stream);
  GtNodeStream *ienv_stream = gt_visitor_stream_new(icnv_stream, ienv);
  gt_queue_add(stream->streams, ienv_stream);
  GtNodeStream *aistream = gt_add_introns_stream_new(ienv_stream);
  gt_queue_add(stream->streams, aistream);
  GtNodeStream *filterstream = agn_filter_stream_new(aistream,typestokeep,NULL);
  gt_queue_add(stream->streams, filterstream);
  stream->in_stream = filterstream;

  return ns;
}

static const GtNodeStreamClass *canon_gene_stream_class(void)
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

static int canon_gene_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *error)
{
  AgnCanonGeneStream *stream = canon_gene_stream_cast(ns);
  return gt_node_stream_next(stream->in_stream, gn, error);
}

static void canon_gene_stream_free(GtNodeStream *ns)
{
  AgnCanonGeneStream *stream = canon_gene_stream_cast(ns);
  while(gt_queue_size(stream->streams) > 0)
  {
    GtNodeStream *s = gt_queue_get(stream->streams);
    gt_node_stream_delete(s);
  }
  gt_queue_delete(stream->streams);
}
