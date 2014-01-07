#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnLocus.h"
#include "AgnLocusStream.h"
#include "AgnTranscriptStream.h"
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
 * @function Feeds feature nodes of type ``locus`` to the output stream.
 */
static int locus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                             GtError *error);

/**
 * @function Determine loci without consideration for the source of each
 * transcript annotation.
 */
static void locus_stream_parse(AgnLocusStream *stream);

/**
 * @function Determine loci while keeping track of which transcripts belong to
 * the reference annotation and which belong to the prediciton annotation.
 */
static void locus_stream_parse_pairwise(AgnLocusStream *stream);

/**
 * @function Given a locus, check the feature index to see if any additional
 * transcripts overlap with this locus.
 */
static int locus_stream_query_overlap(AgnLocusStream *stream, AgnLocus *locus,
                                      GtHashmap *visited);

/**
 * @function Given a locus, check the feature index to see if any additional
 * transcripts overlap with this locus.
 */
static int locus_stream_query_overlap_pairwise(AgnLocusStream *stream,
                                               AgnLocus *locus,
                                               AgnComparisonSource source,
                                               GtHashmap *visited);

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
  stream->in_stream = gt_node_stream_ref(in_stream);
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
  locus_stream_parse(stream);
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
  locus_stream_parse_pairwise(stream);
  stream->out_stream = gt_feature_out_stream_new(stream->loci);

  return ns;
}

bool agn_locus_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  locus_stream_test_data(queue);

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
      GtUword refrtrans = agn_locus_num_refr_transcripts(locus);
      GtUword predtrans = agn_locus_num_pred_transcripts(locus);
      grapetest = grapetest && starts[i] == range.start &&
                  ends[i] == range.end && numrefr[i] == refrtrans &&
                  numpred[i] == predtrans;
      i++;
      //agn_locus_delete(locus);
    }
  }
  agn_unit_test_result(test, "grape test (pairwise)", grapetest);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
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
  gt_node_stream_delete(stream->in_stream);
  gt_node_stream_delete(stream->out_stream);
  if(stream->transcripts != NULL)
    gt_feature_index_delete(stream->transcripts);
  if(stream->refrtrans != NULL)
    gt_feature_index_delete(stream->refrtrans);
  if(stream->predtrans != NULL)
    gt_feature_index_delete(stream->predtrans);
  gt_feature_index_delete(stream->loci);
}

static void locus_stream_parse(AgnLocusStream *stream)
{
  GtUword numseqs, i, j;
  GtError *error = gt_error_new();

  GtStrArray *seqids = gt_feature_index_get_seqids(stream->transcripts, error);
  if(gt_error_is_set(error))
  {
    gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse] error "
                  "retrieving sequence IDs: %s\n", gt_error_get(error));
  }
  numseqs = gt_str_array_size(seqids);

  for(i = 0; i < numseqs; i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtStr *seqidstr = gt_str_new_cstr(seqid);
    GtArray *features;
    features = gt_feature_index_get_features_for_seqid(stream->transcripts,
                                                       seqid, error);
    if(gt_error_is_set(error))
    {
      gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse] "
                    "error retrieving features for sequence '%s': %s", seqid,
                    gt_error_get(error));
    }

    GtHashmap *visited = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    for(j = 0; j < gt_array_size(features); j++)
    {
      GtFeatureNode **transcript = gt_array_get(features, j);
      if(gt_hashmap_get(visited, *transcript) != NULL)
        continue; // Already been processed and assigned to a locus

      gt_hashmap_add(visited, *transcript, *transcript);
      AgnLocus *locus = agn_locus_new(seqidstr);
      agn_locus_add_transcript(locus, *transcript);

      int new_trans_count = 0;
      do
      {
        new_trans_count = locus_stream_query_overlap(stream, locus, visited);
      } while(new_trans_count > 0);
      GtFeatureNode *locusfn = gt_feature_node_cast(locus);
      gt_feature_index_add_feature_node(stream->loci, locusfn, error);
      if(gt_error_is_set(error))
      {
        gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse] "
                      "error adding locus %s[%lu, %lu] to feature index: %s",
                      seqid, gt_genome_node_get_start(locus),
                      gt_genome_node_get_end(locus), gt_error_get(error));
      }
    } // end iterate over features

    gt_hashmap_delete(visited);
    gt_str_delete(seqidstr);
    gt_array_delete(features);
  } // end iterate over seqids

  gt_error_delete(error);
}

static void locus_stream_parse_pairwise(AgnLocusStream *stream)
{
  GtUword numseqs, i, j;
  GtError *error = gt_error_new();

  GtStrArray *refrseqids = gt_feature_index_get_seqids(stream->refrtrans,error);
  GtStrArray *predseqids = gt_feature_index_get_seqids(stream->predtrans,error);
  if(gt_error_is_set(error))
  {
    gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse] error "
                  "retrieving sequence IDs: %s\n", gt_error_get(error));
  }
  GtStrArray *seqids = agn_str_array_union(refrseqids, predseqids);
  numseqs = gt_str_array_size(seqids);

  for(i = 0; i < numseqs; i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtStr *seqidstr = gt_str_new_cstr(seqid);
    GtArray *refr_list;
    refr_list = gt_feature_index_get_features_for_seqid(stream->refrtrans,
                                                        seqid, error);
    if(gt_error_is_set(error))
    {
      gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse_"
                    "pairwise] error retrieving reference features for sequence"
                    " '%s': %s", seqid, gt_error_get(error));
    }

    GtHashmap *visited = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    for(j = 0; j < gt_array_size(refr_list); j++)
    {
      GtFeatureNode **transcript = gt_array_get(refr_list, j);
      if(gt_hashmap_get(visited, *transcript) != NULL)
        continue; // Already been processed and assigned to a locus

      gt_hashmap_add(visited, *transcript, *transcript);
      AgnLocus *locus = agn_locus_new(seqidstr);
      agn_locus_add_refr_transcript(locus, *transcript);

      int new_trans_count = 0;
      do
      {
        int new_refr_count, new_pred_count;
        new_refr_count = locus_stream_query_overlap_pairwise(stream, locus,
                                                             REFERENCESOURCE,
                                                             visited);
        new_pred_count = locus_stream_query_overlap_pairwise(stream, locus,
                                                             PREDICTIONSOURCE,
                                                             visited);
        new_trans_count = new_refr_count + new_pred_count;
      } while(new_trans_count > 0);
      GtFeatureNode *locusfn = gt_feature_node_cast(locus);
      gt_feature_index_add_feature_node(stream->loci, locusfn, error);
      if(gt_error_is_set(error))
      {
        gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse"
                      "_pairwise] error adding locus %s[%lu, %lu] to feature "
                      "index: %s", seqid, gt_genome_node_get_start(locus),
                      gt_genome_node_get_end(locus), gt_error_get(error));
      }
    } // end iterate over features
    gt_array_delete(refr_list);

    GtArray *pred_list;
    pred_list = gt_feature_index_get_features_for_seqid(stream->predtrans,
                                                        seqid, error);
    if(gt_error_is_set(error))
    {
      gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse_"
                    "pairwise] error retrieving prediction features for "
                    "sequence '%s': %s", seqid, gt_error_get(error));
    }
    for(j = 0; j < gt_array_size(pred_list); j++)
    {
      GtFeatureNode **transcript = gt_array_get(pred_list, j);
      if(gt_hashmap_get(visited, *transcript) != NULL)
        continue; // Already been processed and assigned to a locus

      gt_hashmap_add(visited, *transcript, *transcript);
      AgnLocus *locus = agn_locus_new(seqidstr);
      agn_locus_add_pred_transcript(locus, *transcript);

      int new_trans_count = 0;
      do
      {
        new_trans_count = locus_stream_query_overlap_pairwise(stream, locus,
                                                              PREDICTIONSOURCE,
                                                              visited);
      } while(new_trans_count > 0);
      GtFeatureNode *locusfn = gt_feature_node_cast(locus);
      gt_feature_index_add_feature_node(stream->loci, locusfn, error);
      if(gt_error_is_set(error))
      {
        gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_parse"
                      "_pairwise] error adding locus %s[%lu, %lu] to feature "
                      "index: %s", seqid, gt_genome_node_get_start(locus),
                      gt_genome_node_get_end(locus), gt_error_get(error));
      }
    } // end iterate over features
    gt_array_delete(pred_list);
    gt_str_delete(seqidstr);
    gt_hashmap_delete(visited);

  } // end iterate over seqids

  gt_str_array_delete(seqids);
  gt_error_delete(error);
}

static int locus_stream_query_overlap(AgnLocusStream *stream, AgnLocus *locus,
                                      GtHashmap *visited)
{
  GtError *error = gt_error_new();
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange range = gt_genome_node_get_range(locus);
  bool has_seqid;
  gt_feature_index_has_seqid(stream->transcripts, &has_seqid, gt_str_get(seqid),
                             error);
  gt_assert(has_seqid);

  int new_trans_count = 0;
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  gt_feature_index_get_features_for_range(stream->transcripts, overlapping,
                                          gt_str_get(seqid), &range, error);
  if(gt_error_is_set(error))
  {
    gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_query_overlap]"
                  "error retrieving overlapping transcripts for locus %s[%lu, "
                  "%lu]: %s\n", gt_str_get(seqid), range.start, range.end,
                  gt_error_get(error));
  }
  gt_error_delete(error);

  while(gt_array_size(overlapping) > 0)
  {
    GtFeatureNode **fn = gt_array_pop(overlapping);
    if(gt_hashmap_get(visited, *fn) == NULL)
    {
      gt_hashmap_add(visited, *fn, *fn);
      agn_locus_add_transcript(locus, *fn);
      new_trans_count++;
    }
  }
  gt_array_delete(overlapping);

  return new_trans_count;
}

static int locus_stream_query_overlap_pairwise(AgnLocusStream *stream,
                                               AgnLocus *locus,
                                               AgnComparisonSource source,
                                               GtHashmap *visited)
{
  GtError *error = gt_error_new();
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange range = gt_genome_node_get_range(locus);
  bool has_seqid;
  GtFeatureIndex *transcripts = (source == REFERENCESOURCE) ? stream->refrtrans
                                                            : stream->predtrans;
  gt_feature_index_has_seqid(transcripts, &has_seqid, gt_str_get(seqid), error);
  if(!has_seqid)
  {
    gt_error_delete(error);
    return 0;
  }

  int new_trans_count = 0;
  GtArray *overlapping = gt_array_new( sizeof(GtFeatureNode *) );
  gt_feature_index_get_features_for_range(transcripts, overlapping,
                                          gt_str_get(seqid), &range, error);
  if(gt_error_is_set(error))
  {
    const char *src = (source == REFERENCESOURCE) ? "reference" : "prediction";
    gt_logger_log(stream->logger, "[AgnLocusStream::locus_stream_query_overlap"
                  "_pairwise] error retrieving overlapping %s transcripts for "
                  "locus %s[%lu, %lu]: %s\n", src, gt_str_get(seqid),
                  range.start, range.end, gt_error_get(error));
  }

  while(gt_array_size(overlapping) > 0)
  {
    GtFeatureNode **fn = gt_array_pop(overlapping);
    if(gt_hashmap_get(visited, *fn) == NULL)
    {
      gt_hashmap_add(visited, *fn, *fn);
      agn_locus_add(locus, *fn, source);
      new_trans_count++;
    }
  }

  gt_array_delete(overlapping);
  gt_error_delete(error);
  return new_trans_count;
}

static void locus_stream_test_data(GtQueue *queue)
{
  GtError *error = gt_error_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);

  const char *file = "data/gff3/grape-refr-mrnas.gff3";
  GtNodeStream *refrgff3 = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)refrgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)refrgff3);
  GtNodeStream *refrtrans = agn_transcript_stream_new(refrgff3, logger);

  file = "data/gff3/grape-pred-mrnas.gff3";
  GtNodeStream *predgff3 = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)predgff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)predgff3);
  GtNodeStream *predtrans = agn_transcript_stream_new(predgff3, logger);

  GtNodeStream *locusstream;
  locusstream = agn_locus_stream_new_pairwise(refrtrans, predtrans, logger);
  GtArray *loci = gt_array_new( sizeof(AgnLocus *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(locusstream, loci, error);
  gt_node_stream_pull(arraystream, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[AgnLocusStream::locus_stream_test_data] error processing "
            "node stream: %s\n", gt_error_get(error));
  }
  gt_array_reverse(loci);
  while(gt_array_size(loci) > 0)
  {
    AgnLocus **locus = gt_array_pop(loci);
    gt_queue_add(queue, *locus);
  }
  gt_array_delete(loci);

  gt_node_stream_delete(refrgff3);
  gt_node_stream_delete(predgff3);
  gt_node_stream_delete(refrtrans);
  gt_node_stream_delete(predtrans);
  gt_node_stream_delete(locusstream);
  //gt_node_stream_delete(arraystream);
  gt_error_delete(error);
  gt_logger_delete(logger);
}

