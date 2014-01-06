#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnLocusStream.h"
#include "AgnTypecheck.h"

#define locus_stream_cast(GS)\
        gt_node_stream_cast(locus_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnLocusStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtFeatureIndex *transcripts;
  GtFeatureIndex *refrtrans;
  GtFeatureIndex *predtrans;
  GtFeatureIndex *loci;
  GtNodeStream *out_stream;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Function that implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* locus_stream_class(void);

/**
 * @function Class destructor.
 */
static void locus_stream_free(GtNodeStream *ns);

/**
 * @function Determine loci without consideration for the source of each
 * transcript annotation.
 */
static void
locus_stream_parse(AgnLocusStream *stream, GtLogger *logger);

/**
 * @function Determine loci while keeping track of which transcripts belong to
 * the reference annotation and which belong to the prediciton annotation.
 */
static void
locus_stream_parse_pairwise(AgnLocusStream *stream, GtLogger *logger);

/**
 * @function Feeds feature nodes of type ``locus`` to the output stream.
 */
static int locus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                             GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void locus_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_locus_stream_new(GtNodeStream *in_stream, GtLogger *logger)
{
  gt_assert(in_stream);

  GtNodeStream *ns = gt_node_stream_create(locus_stream_class(), false);
  AgnLocusStream *stream = locus_stream_cast(ns);
  stream->logger = logger;

  GtError *error = gt_error_new();
  stream->transcripts = gt_feature_index_memory_new();
  stream->refrtrans = NULL;
  stream->predtrans = NULL;
  GtNodeStream *trans_stream = gt_feature_in_stream_new(in_stream,
                                                        stream->transcripts);
  int result = gt_node_stream_pull(trans_stream, error);
  if(result == -1)
  {
    gt_assert(gt_error_is_set(error));
    gt_logger_log(logger, "[AgnLocusStream::agn_locus_stream_new] error "
                  "processing input: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(trans_stream);
  gt_error_delete(error);
  
  stream->loci = gt_feature_index_memory_new();
  locus_stream_parse(stream, logger);
  stream->out_stream = gt_feature_out_stream_new(stream->loci);

  return ns;
}

GtNodeStream *agn_locus_stream_new_pairwise(GtNodeStream *refr_stream,
                                            GtNodeStream *pred_stream,
                                            GtLogger *logger)
{
  gt_assert(refr_stream && pred_stream);

  GtNodeStream *ns = gt_node_stream_create(locus_stream_class(), false);
  AgnLocusStream *stream = locus_stream_cast(ns);
  stream->logger = logger;

  GtError *error = gt_error_new();
  stream->transcripts = NULL;
  stream->refrtrans = gt_feature_index_memory_new();
  stream->predtrans = gt_feature_index_memory_new();
  GtNodeStream *refr_instream = gt_feature_in_stream_new(refr_stream,
                                                         stream->refrtrans);
  int result = gt_node_stream_pull(refr_instream, error);
  if(result == -1)
  {
    gt_assert(gt_error_is_set(error));
    gt_logger_log(logger, "[AgnLocusStream::agn_locus_stream_new_pairwise] "
                 "error processing reference input: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(refr_instream);
  gt_error_delete(error);
  GtNodeStream *pred_instream = gt_feature_in_stream_new(pred_stream,
                                                         stream->predtrans);
  result = gt_node_stream_pull(pred_instream, error);
  if(result == -1)
  {
    gt_assert(gt_error_is_set(error));
    gt_logger_log(logger, "[AgnLocusStream::agn_locus_stream_new_pairwise] "
                 "error processing prediction input: %s\n",gt_error_get(error));
  }
  gt_node_stream_delete(pred_instream);
  gt_error_delete(error);
  
  stream->loci = gt_feature_index_memory_new();
  locus_stream_parse_pairwise(stream, logger);
  stream->out_stream = gt_feature_out_stream_new(stream->loci);

  return ns;
}

bool agn_locus_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  locus_stream_test_data(queue);

  gt_queue_delete(queue);
  // return agn_unit_test_success(test);
  return false;
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

static int locus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                             GtError *error)
{
  AgnLocusStream *stream = stream = locus_stream_cast(ns);
  return gt_node_stream_next(stream->out_stream, gn, error);
}

static void locus_stream_free(GtNodeStream *ns)
{
  AgnLocusStream *stream = locus_stream_cast(ns);
  gt_node_stream_delete(stream->out_stream);
  if(stream->transcripts != NULL)
    gt_feature_index_delete(stream->transcripts);
  if(stream->refrtrans != NULL)
    gt_feature_index_delete(stream->refrtrans);
  if(stream->predtrans != NULL)
    gt_feature_index_delete(stream->predtrans);
  gt_feature_index_delete(stream->loci);
}

static void
locus_stream_parse(AgnLocusStream *stream, GtLogger *logger)
{
  
}

static void
locus_stream_parse_pairwise(AgnLocusStream *stream, GtLogger *logger)
{
  
}

static void locus_stream_test_data(GtQueue *queue)
{
  
}