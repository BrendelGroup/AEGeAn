#include "core/hashmap_api.h"
#include "core/logger_api.h"
#include "core/queue_api.h"
#include "AgnFilterStream.h"
#include "AgnGeneStream.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
#include "AgnUtils.h"
#include "AgnTypecheck.h"

#define gene_stream_cast(GS)\
        gt_node_stream_cast(gene_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnGeneStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtQueue *streams;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Function that implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* gene_stream_class(void);

/**
 * @function Class destructor.
 */
static void gene_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they pass validation.
 */
static int gene_stream_next(GtNodeStream *ns, GtGenomeNode **gn,GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void gene_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------


GtNodeStream* agn_gene_stream_new(GtNodeStream *in_stream, GtLogger *logger)
{
  GtNodeStream *ns;
  AgnGeneStream *stream;
  gt_assert(in_stream);

  ns = gt_node_stream_create(gene_stream_class(), false);
  stream = gene_stream_cast(ns);
  stream->logger = logger;
  stream->streams = gt_queue_new();
  gt_queue_add(stream->streams, gt_node_stream_ref(in_stream));

  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                          gt_free_func);
  gt_hashmap_add(typestokeep, gt_cstr_dup("gene"), gt_cstr_dup("gene"));

  GtNodeStream *ic_stream = agn_infer_cds_stream_new(in_stream, logger);
  gt_queue_add(stream->streams, ic_stream);
  GtNodeStream *ie_stream = agn_infer_exons_stream_new(ic_stream, logger);
  gt_queue_add(stream->streams, ie_stream);
  GtNodeStream *filterstream = agn_filter_stream_new(ie_stream, typestokeep);
  gt_queue_add(stream->streams, filterstream);
  stream->in_stream = filterstream;

  gt_hashmap_delete(typestokeep);
  return ns;
}

bool agn_gene_stream_unit_test(AgnUnitTest *test)
{
  gene_stream_test_data(NULL);
  return false;
}

static const GtNodeStreamClass *gene_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnGeneStream),
                                   gene_stream_free,
                                   gene_stream_next);
  }
  return nsc;
}

static void gene_stream_free(GtNodeStream *ns)
{
  AgnGeneStream *stream = gene_stream_cast(ns);
  while(gt_queue_size(stream->streams) > 0)
  {
    GtNodeStream *s = gt_queue_get(stream->streams);
    gt_node_stream_delete(s);
  }
  gt_queue_delete(stream->streams);
}

static int gene_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *error)
{
  AgnGeneStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = gene_stream_cast(ns);

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

    gt_assert(agn_typecheck_gene(fn));

    GtUword num_valid_mrnas = 0;
    GtFeatureNode *current;
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    GtQueue *invalid_mrnas = gt_queue_new();
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      if(!agn_typecheck_mrna(current))
        continue;
      
      GtArray *cds     = agn_typecheck_select(current, agn_typecheck_cds);
      GtArray *exons   = agn_typecheck_select(current, agn_typecheck_exon);
      GtArray *introns = agn_typecheck_select(current, agn_typecheck_intron);

      bool keepmrna = true;
      if(gt_array_size(cds) < 1)
      {
        const char *mrnaid = gt_feature_node_get_attribute(current, "ID");
        gt_logger_log(stream->logger, "ignoring mRNA '%s': no CDS", mrnaid);
        keepmrna = false;
      }
      if(gt_array_size(exons) != gt_array_size(introns) + 1)
      {
        const char *mrnaid = gt_feature_node_get_attribute(current, "ID");
        gt_logger_log(stream->logger, "error: mRNA '%s' has %lu exons but "
                      "%lu introns", mrnaid, gt_array_size(exons),
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
    gt_feature_node_iterator_delete(iter);
    while(gt_queue_size(invalid_mrnas) > 0)
    {
      GtFeatureNode *mrna = gt_queue_get(invalid_mrnas);
      agn_feature_node_remove_tree(fn, mrna);
    }
    gt_queue_delete(invalid_mrnas);

    if(num_valid_mrnas > 0)
      return 0;
    else
    {
      char idstr[1024];
      const char *id = gt_feature_node_get_attribute(fn, "ID");
      if(id == NULL)
      {
        GtStr *seqid = gt_genome_node_get_seqid(*gn);
        GtRange rng = gt_genome_node_get_range(*gn);
        sprintf(idstr, "%s[%lu, %lu]", gt_str_get(seqid), rng.start, rng.end);
        id = idstr;
      }
      gt_logger_log(stream->logger, "warning: found no valid mRNAs for gene "
                    "'%s'", id);
      gt_genome_node_delete(*gn);
    }
  }

  return 0;
}

static void gene_stream_test_data(GtQueue *queue)
{
  
}
