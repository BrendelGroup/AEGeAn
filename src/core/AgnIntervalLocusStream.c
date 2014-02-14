#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnIntervalLocusStream.h"
#include "AgnLocus.h"
#include "AgnLocusStream.h"
#include "AgnTranscriptStream.h"

#define ilocus_stream_cast(GS)\
        gt_node_stream_cast(ilocus_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnIntervalLocusStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *ilocus_stream;
  GtNodeStream *locus_stream;
  GtNodeStream *fstream;
  GtFeatureIndex *in_loci;
  GtFeatureIndex *out_loci;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Function that implements the GtNodeStream interface for this class.
 */
static const GtNodeStreamClass* ilocus_stream_class(void);

/**
 * @function Class destructor.
 */
static void ilocus_stream_free(GtNodeStream *ns);

/**
 * @function Feeds feature nodes of type ``locus`` to the output stream.
 */
static int ilocus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error);

/**
 * @function Use gene/transcript loci to determine iLoci using the given delta.
 */
static void ilocus_stream_parse(AgnIntervalLocusStream *stream, GtUword delta,
                                bool skipterminal);

/**
 * @function Handle the first locus while parsing iLoci.
 */
static void ilocus_stream_parse_initial(AgnIntervalLocusStream *stream,
                                        GtUword delta, bool skipterminal,
                                        AgnSequenceRegion *region,
                                        GtArray *seqloci, GtError *error);

/**
 * @function Handle internal loci while parsing iLoci.
 */
static void ilocus_stream_parse_internal(AgnIntervalLocusStream *stream,
                                         GtUword delta, bool skipterminal,
                                         GtStr *seqidstr, GtArray *seqloci,
                                         GtError *error);

/**
 * @function Parse all iLoci for a given sequence.
 */
static void ilocus_stream_parse_seq(AgnIntervalLocusStream *stream,
                                    GtUword delta, bool skipterminal,
                                    const char *seqid, GtError *error);

/**
 * @function Handle the final locus while parsing iLoci.
 */
static void ilocus_stream_parse_terminal(AgnIntervalLocusStream *stream,
                                         GtUword delta, bool skipterminal,
                                         AgnSequenceRegion *region,
                                         GtArray *seqloci, GtError *error);

/**
 * @function Generate data for unit testing.
 */
static GtNodeStream*
ilocus_stream_test_data(GtFeatureIndex *iloci,GtNodeStream *s1,GtNodeStream *s2,
                        GtUword delta, bool skipends, GtLogger *logger);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream *agn_interval_locus_stream_new(GtNodeStream *locus_stream,
                                            GtUword delta,
                                            bool skipterminal,
                                            GtLogger *logger)
{
  gt_assert(locus_stream);

  GtNodeStream *ns = gt_node_stream_create(ilocus_stream_class(), false);
  AgnIntervalLocusStream *stream = ilocus_stream_cast(ns);
  stream->logger = logger;

  GtError *error = gt_error_new();
  stream->locus_stream = gt_node_stream_ref(locus_stream);
  stream->in_loci = gt_feature_index_memory_new();
  stream->fstream = gt_feature_out_stream_new(locus_stream, stream->in_loci);
  int result = gt_node_stream_pull(stream->fstream, error);
  if(result == -1)
  {
    gt_assert(gt_error_is_set(error));
    gt_logger_log(logger, "[AgnIntervalLocusStream::agn_interval_locus_stream_"
                  "new] error processing input: %s\n", gt_error_get(error));
  }

  stream->out_loci = gt_feature_index_memory_new();
  agn_feature_index_copy_regions(stream->out_loci, stream->in_loci, true,error);
  gt_error_delete(error);

  ilocus_stream_parse(stream, delta, skipterminal);
  stream->ilocus_stream = gt_feature_in_stream_new(stream->out_loci);
  gt_feature_in_stream_use_orig_ranges((GtFeatureInStream *)stream->out_loci);

  return ns;
}

bool agn_interval_locus_stream_unit_test(AgnUnitTest *test)
{
  GtFeatureIndex *iloci = gt_feature_index_memory_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  const char *infile = "data/gff3/ilocus.in.gff3";
  GtNodeStream *gff3 = gt_gff3_in_stream_new_unsorted(1, &infile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3);
  GtNodeStream *fstream = ilocus_stream_test_data(iloci, gff3, NULL, 200, false,
                                                  logger);

  GtStrArray *seqids = gt_feature_index_get_seqids(iloci, error);
  if(gt_str_array_size(seqids) != 25)
  {
    agn_unit_test_result(test, "Test 1", false);
    return agn_unit_test_success(test);
  }

  const char *seqid = gt_str_array_get(seqids, 0);
  GtArray *seqloci = gt_feature_index_get_features_for_seqid(iloci,seqid,error);
  bool test1_1 = (gt_array_size(seqloci) == 1);
  if(test1_1)
  {
    AgnLocus *locus = *(AgnLocus **)gt_array_get(seqloci, 0);
    GtRange range = gt_genome_node_get_range(locus);
    test1_1 = range.start == 1 && range.end == 900;
  }
  agn_unit_test_result(test, "Test 1.1", test1_1);
  gt_array_delete(seqloci);

  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    seqloci = gt_feature_index_get_features_for_seqid(iloci, seqid, error);
    while(gt_array_size(seqloci) > 0)
    {
      AgnLocus **locus = gt_array_pop(seqloci);
      gt_genome_node_delete(*locus);
    }
    gt_array_delete(seqloci);
  }

  gt_feature_index_delete(iloci);
  gt_str_array_delete(seqids);
  gt_logger_delete(logger);
  gt_node_stream_delete(gff3);
  gt_node_stream_delete(fstream);
  gt_error_delete(error);
  return agn_unit_test_success(test);
}

static const GtNodeStreamClass *ilocus_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if(!nsc)
  {
    nsc = gt_node_stream_class_new(sizeof (AgnIntervalLocusStream),
                                   ilocus_stream_free,
                                   ilocus_stream_next);
  }
  return nsc;
}

static void ilocus_stream_free(GtNodeStream *ns)
{
  AgnIntervalLocusStream *stream = ilocus_stream_cast(ns);
  gt_node_stream_delete(stream->ilocus_stream);
  gt_feature_index_delete(stream->in_loci);
  gt_feature_index_delete(stream->out_loci);
  gt_node_stream_delete(stream->locus_stream);
  gt_node_stream_delete(stream->fstream);
}

static int ilocus_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                              GtError *error)
{
  AgnIntervalLocusStream *stream = ilocus_stream_cast(ns);
  return gt_node_stream_next(stream->ilocus_stream, gn, error);
}

static void ilocus_stream_parse(AgnIntervalLocusStream *stream, GtUword delta,
                                bool skipterminal)
{
  GtError *error = gt_error_new();
  GtStrArray *seqids = gt_feature_index_get_seqids(stream->in_loci, error);
  if(gt_error_is_set(error))
  {
    gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_stream_"
                  "parse] error retrieving sequence IDs: %s\n",
                  gt_error_get(error));
  }

  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    ilocus_stream_parse_seq(stream, delta, skipterminal, seqid, error);
  }

  gt_error_delete(error);
  gt_str_array_delete(seqids);
}

static void ilocus_stream_parse_initial(AgnIntervalLocusStream *stream,
                                        GtUword delta, bool skipterminal,
                                        AgnSequenceRegion *region,
                                        GtArray *seqloci, GtError *error)
{
  GtRange *seqrange = &region->range;
  AgnLocus *l1 = *(AgnLocus **)gt_array_get(seqloci, 0);
  GtRange r1 = gt_genome_node_get_range(l1);
  if(r1.start >= seqrange->start + (2*delta))
  {
    if(gt_array_size(seqloci) == 1)
      agn_locus_set_range(l1, r1.start - delta, r1.end);

    if(!skipterminal)
    {
      AgnLocus *ilocus = agn_locus_new(region->seqid);
      agn_locus_set_range(ilocus, seqrange->start, r1.start - delta - 1);
      GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
      gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
    }
  }
  else
  {
    agn_locus_set_range(l1, seqrange->start, r1.end);
  }
}

static void ilocus_stream_parse_internal(AgnIntervalLocusStream *stream,
                                         GtUword delta, bool skipterminal,
                                         GtStr *seqidstr, GtArray *seqloci,
                                         GtError *error)
{
  GtUword i;
  for(i = 0; i < gt_array_size(seqloci) - 1; i++)
  {
    AgnLocus *llocus = *(AgnLocus **)gt_array_get(seqloci, i);
    AgnLocus *rlocus = *(AgnLocus **)gt_array_get(seqloci, i+1);
    GtRange lrange = gt_genome_node_get_range(llocus);
    GtRange rrange = gt_genome_node_get_range(rlocus);

    GtFeatureNode *llocusfn = gt_feature_node_cast(llocus);
    gt_feature_index_add_feature_node(stream->out_loci, llocusfn, error);
    if(lrange.end + delta >= rrange.start)
    {
      agn_locus_set_range(llocus, lrange.start, rrange.start - 1);
      agn_locus_set_range(rlocus, lrange.end + 1, rrange.end);
    }
    else if(lrange.end + (2*delta) >= rrange.start)
    {
      agn_locus_set_range(llocus, lrange.start, lrange.end + delta);
      agn_locus_set_range(rlocus, rrange.start - delta, rrange.end);
    }
    else if(lrange.end + (3*delta) >= rrange.start)
    {
      GtUword midpoint = (lrange.end + rrange.start) / 2;
      agn_locus_set_range(llocus, lrange.start, midpoint);
      agn_locus_set_range(rlocus, midpoint + 1, rrange.end);
    }
    else
    {
      agn_locus_set_range(llocus, lrange.start, lrange.end + delta);
      agn_locus_set_range(rlocus, rrange.start - delta, rrange.end);

      AgnLocus *ilocus = agn_locus_new(seqidstr);
      agn_locus_set_range(ilocus, lrange.end + delta + 1,
                          rrange.start - delta - 1);
      GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
      gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
    }

    if(i == gt_array_size(seqloci) - 2)
    {
      GtFeatureNode *rlocusfn = gt_feature_node_cast(rlocus);
      gt_feature_index_add_feature_node(stream->out_loci, rlocusfn, error);
    }
  }
}

static void ilocus_stream_parse_seq(AgnIntervalLocusStream *stream,
                                    GtUword delta, bool skipterminal,
                                    const char *seqid, GtError *error)
{
  GtStr *seqidstr = gt_str_new_cstr(seqid);
  GtRange seqrange;
  gt_feature_index_get_orig_range_for_seqid(stream->in_loci, &seqrange, seqid,
                                            error);
  AgnSequenceRegion seqregion = { seqidstr, seqrange };
  GtArray *seqloci = gt_feature_index_get_features_for_seqid(stream->in_loci,
                                                             seqid, error);
  if(seqloci == NULL)
    return;
  GtUword nloci = gt_array_size(seqloci);
  if(nloci > 1)
    gt_array_sort(seqloci, (GtCompare)agn_genome_node_compare);

  // Handle trivial case
  if(nloci == 0)
  {
    AgnLocus *ilocus = agn_locus_new(seqidstr);
    GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
    agn_locus_set_range(ilocus, seqrange.start, seqrange.end);
    gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
    gt_array_delete(seqloci);
    gt_str_delete(seqidstr);
    return;
  }

  ilocus_stream_parse_initial(stream, delta, skipterminal,
                              &seqregion, seqloci, error);

  ilocus_stream_parse_internal(stream, delta, skipterminal,
                               seqidstr, seqloci, error);

  ilocus_stream_parse_terminal(stream, delta, skipterminal,
                               &seqregion, seqloci, error);

  gt_str_delete(seqidstr);
  gt_array_delete(seqloci);
}

static void ilocus_stream_parse_terminal(AgnIntervalLocusStream *stream,
                                         GtUword delta, bool skipterminal,
                                         AgnSequenceRegion *region,
                                         GtArray *seqloci, GtError *error)
{
  GtRange *seqrange = &region->range;
  GtUword nloci = gt_array_size(seqloci);
  AgnLocus *locus = *(AgnLocus **)gt_array_get(seqloci, nloci - 1);
  GtRange range = gt_genome_node_get_range(locus);
  if(range.end <= seqrange->end - (2*delta))
  {
    agn_locus_set_range(locus, range.start, range.end + delta);
  
    if(!skipterminal)
    {
      AgnLocus *ilocus = agn_locus_new(region->seqid);
      agn_locus_set_range(ilocus, range.end + delta + 1, seqrange->end);
      GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
      gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
    }
  }
  else
  {
    agn_locus_set_range(locus, range.start, seqrange->end);
  }

  if(nloci == 1)
  {
    GtFeatureNode *locusfn = gt_feature_node_cast(locus);
    gt_feature_index_add_feature_node(stream->out_loci, locusfn, error);
  }
}

static GtNodeStream*
ilocus_stream_test_data(GtFeatureIndex *iloci,GtNodeStream *s1,GtNodeStream *s2,
                        GtUword delta, bool skipends, GtLogger *logger)
{
  gt_assert(s1);
  GtError *error = gt_error_new();

  GtNodeStream *locusstream, *ilocusstream;
  if(s2 == NULL)
    locusstream = agn_locus_stream_new(s1, logger);
  else
    locusstream = agn_locus_stream_new_pairwise(s1, s2, logger);
  ilocusstream = agn_interval_locus_stream_new(locusstream, delta, skipends,
                                               logger);

  GtNodeStream *fstream = gt_feature_out_stream_new(ilocusstream, iloci);
  gt_node_stream_pull(fstream, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[AgnLocusStream::interval_locus_stream_test_data] error "
            "processing node stream: %s\n", gt_error_get(error));
  }

  gt_node_stream_delete(locusstream);
  gt_node_stream_delete(ilocusstream);
  gt_error_delete(error);
  return fstream;
}
