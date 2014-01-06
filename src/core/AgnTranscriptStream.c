#include "core/array_api.h"
#include "core/hashmap_api.h"
#include "core/queue_api.h"
#include "AgnFilterStream.h"
#include "AgnInferCDSVisitor.h"
#include "AgnInferExonsVisitor.h"
#include "AgnTranscriptStream.h"
#include "AgnUtils.h"
#include "AgnTypecheck.h"

#define transcript_stream_cast(GS)\
        gt_node_stream_cast(transcript_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnTranscriptStream
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
static const GtNodeStreamClass* transcript_stream_class(void);

/**
 * @function Class destructor.
 */
static void transcript_stream_free(GtNodeStream *ns);

/**
 * @function Pulls nodes from the input stream and feeds them to the output
 * stream if they pass validation.
 */
static int transcript_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *error);

/**
 * @function Generate data for unit testing.
 */
static void transcript_stream_test_data(GtQueue *queue);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_transcript_stream_new(GtNodeStream *in_stream,
                                        GtLogger *logger)
{
  GtNodeStream *ns;
  AgnTranscriptStream *stream;
  gt_assert(in_stream);

  ns = gt_node_stream_create(transcript_stream_class(), false);
  stream = transcript_stream_cast(ns);
  stream->logger = logger;
  stream->streams = gt_queue_new();
  gt_queue_add(stream->streams, gt_node_stream_ref(in_stream));

  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  gt_hashmap_add(typestokeep, "mRNA", "mRNA");
  gt_hashmap_add(typestokeep, "rRNA", "rRNA");
  gt_hashmap_add(typestokeep, "tRNA", "tRNA");

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

bool agn_transcript_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  transcript_stream_test_data(queue);
  gt_assert(gt_queue_size(queue) == 4);

  GtFeatureNode *fn = gt_queue_get(queue);
  GtArray *utr3p = agn_typecheck_select(fn, agn_typecheck_utr3p);
  GtArray *utr5p = agn_typecheck_select(fn, agn_typecheck_utr5p);
  GtRange range = gt_genome_node_get_range((GtGenomeNode *)fn);
  bool test1 = (gt_array_size(utr3p) == 1 && gt_array_size(utr5p) == 0);
  if(test)
  {
    GtGenomeNode *utr = *(GtGenomeNode **)gt_array_get(utr3p, 0);
    GtRange utrrange = gt_genome_node_get_range(utr);
    test1 = (utrrange.start == 2926 && utrrange.end == 3035 &&
             range.start == 600 && range.end == 3035);
  }
  agn_unit_test_result(test, "mRNA only", test1);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(utr3p);
  gt_array_delete(utr5p);

  GtGenomeNode *fn1 = gt_queue_get(queue);
  GtGenomeNode *fn2 = gt_queue_get(queue);
  GtRange range1 = gt_genome_node_get_range(fn1);
  GtRange range2 = gt_genome_node_get_range(fn2);
  bool test2 = (range1.start == 5000 && range1.end == 9000 &&
                range2.start == 6000 && range2.end == 9000);
  agn_unit_test_result(test, "1 gene, 2 mRNAs", test2);
  gt_genome_node_delete(fn1);
  gt_genome_node_delete(fn2);

  fn = gt_queue_get(queue);
  range = gt_genome_node_get_range((GtGenomeNode *)fn);
  bool test3 = (range.start == 15000 && range.end == 19000);
  agn_unit_test_result(test, "1 gene, 2 mRNAs (one sans CDS)", test3);
  gt_genome_node_delete((GtGenomeNode *)fn);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static const GtNodeStreamClass *transcript_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnTranscriptStream),
                                   transcript_stream_free,
                                   transcript_stream_next);
  }
  return nsc;
}

static int transcript_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *error)
{
  AgnTranscriptStream *stream;
  GtFeatureNode *fn;
  int had_err;
  gt_error_check(error);
  stream = transcript_stream_cast(ns);

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
      if(agn_typecheck_mrna(current))
      {
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
    }
    gt_feature_node_iterator_delete(iter);
    while(gt_queue_size(invalid_mrnas) > 0)
    {
      GtGenomeNode *mrna = gt_queue_get(invalid_mrnas);
      gt_genome_node_delete(mrna);
    }
    gt_queue_delete(invalid_mrnas);

    if(num_valid_mrnas > 0)
      return 0;
  }

  return 0;
}

static void transcript_stream_free(GtNodeStream *ns)
{
  AgnTranscriptStream *stream = transcript_stream_cast(ns);
  while(gt_queue_size(stream->streams) > 0)
  {
    GtNodeStream *s = gt_queue_get(stream->streams);
    gt_node_stream_delete(s);
  }
  gt_queue_delete(stream->streams);
}

static void transcript_stream_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  const char *file = "data/gff3/trans-stream-data.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  FILE *log = fopen("/dev/null", "w");
  if(log == NULL)
  {
    fprintf(stderr, "[AgnTranscriptStream::transcript_stream_test_data] error "
            "opening /dev/null");
  }
  GtLogger *logger = gt_logger_new(true, "", log);
  GtNodeStream *stream = agn_transcript_stream_new(gff3in, logger);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(stream, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnTranscriptStream::transcript_stream_test_data] error "
            "processing features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(stream);
  gt_node_stream_delete(arraystream);
  gt_logger_delete(logger);
  fclose(log);
  gt_array_sort(feats, (GtCompare)agn_genome_node_compare);
  gt_array_reverse(feats);
  while(gt_array_size(feats) > 0)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_pop(feats);
    gt_queue_add(queue, fn);
  }
  gt_array_delete(feats);
  gt_error_delete(error);
}