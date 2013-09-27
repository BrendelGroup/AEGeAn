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
  AgnLogger *logger;
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
  stream->logger = logger;

  GtNodeVisitor *icnv = agn_infer_cds_visitor_new(logger);
  GtNodeVisitor *ienv = agn_infer_exons_visitor_new(logger);
  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "gene", "gene");

  gt_queue_add(stream->streams, gt_node_stream_ref(in_stream));
  GtNodeStream *icnv_stream = gt_visitor_stream_new(in_stream, icnv);
  gt_queue_add(stream->streams, icnv_stream);
  GtNodeStream *ienv_stream = gt_visitor_stream_new(icnv_stream, ienv);
  gt_queue_add(stream->streams, ienv_stream);
  GtNodeStream *filterstream = agn_filter_stream_new(ienv_stream, typestokeep);
  gt_queue_add(stream->streams, filterstream);
  stream->in_stream = filterstream;

  gt_hashmap_delete(typestokeep);
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
  AgnCanonGeneStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = canon_gene_stream_cast(ns);

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

    GtUword num_valid_mrnas = 0;
    GtFeatureNode *current;
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    GtQueue *invalid_mrnas = gt_queue_new();
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      if(agn_gt_feature_node_is_mrna_feature(current))
      {
        GtArray *cds     = agn_gt_feature_node_children_of_type(current,
                                            agn_gt_feature_node_is_cds_feature);
        GtArray *exons   = agn_gt_feature_node_children_of_type(current,
                                           agn_gt_feature_node_is_exon_feature);
        GtArray *introns = agn_gt_feature_node_children_of_type(current,
                                         agn_gt_feature_node_is_intron_feature);

        bool keepmrna = true;
        if(gt_array_size(cds) < 1)
        {
          const char *mrnaid = gt_feature_node_get_attribute(current, "ID");
          agn_logger_log_warning(stream->logger, "ignoring mRNA '%s': no CDS",
                                 mrnaid);
          keepmrna = false;
        }
        if(gt_array_size(exons) != gt_array_size(introns) + 1)
        {
          const char *mrnaid = gt_feature_node_get_attribute(current, "ID");
          agn_logger_log_error(stream->logger, "error: mRNA '%s' has %lu exons "
                               "but %lu introns", mrnaid, gt_array_size(exons),
                               gt_array_size(introns));
          keepmrna = false;
        }

        if(keepmrna)
          num_valid_mrnas++;
        else
          gt_queue_add(invalid_mrnas, current);

        gt_array_delete(cds);
        gt_array_delete(exons);
        gt_array_delete(introns);
      }
    }
    gt_feature_node_iterator_delete(iter);
    gt_queue_delete(invalid_mrnas);

    if(num_valid_mrnas > 0)
      return 0;
  }

  return 0;
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
