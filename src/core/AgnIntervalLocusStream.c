#include "core/queue_api.h"
#include "extended/feature_index_memory_api.h"
#include "AgnLocus.h"
#include "AgnIntervalLocusStream.h"

#define ilocus_stream_cast(GS)\
        gt_node_stream_cast(ilocus_stream_class(), GS)

//------------------------------------------------------------------------------
// Data structure definition
//------------------------------------------------------------------------------
struct AgnIntervalLocusStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *ilocus_stream;
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
  AgnIntervalLocusStream *stream = ilocus_stream_cast(locus_stream);
  stream->logger = logger;

  GtError *error = gt_error_new();
  stream->in_loci = gt_feature_index_memory_new();
  GtNodeStream *fstream = gt_feature_out_stream_new(locus_stream,
                                                    stream->in_loci);
  int result = gt_node_stream_pull(fstream, error);
  if(result == -1)
  {
    gt_assert(gt_error_is_set(error));
    gt_logger_log(logger, "[AgnIntervalLocusStream::agn_interval_locus_stream_"
                  "new] error processing input: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(locus_stream);
  gt_node_stream_delete(fstream);
  gt_error_delete(error);

  stream->out_loci = gt_feature_index_memory_new();
  ilocus_stream_parse(stream, delta, skipterminal);
  stream->ilocus_stream = gt_feature_in_stream_new(stream->out_loci);

  return ns;
}

bool agn_interval_locus_stream_unit_test(AgnUnitTest *test)
{
  return false;
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
    GtStr *seqidstr = gt_str_new_cstr(seqid);
    GtArray *seqloci = gt_feature_index_get_features_for_seqid(stream->in_loci,
                                                               seqid,
                                                               error);
    if(gt_error_is_set(error))
    {
      gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_stream_"
                    "parse] error retrieving features for sequence '%s': %s\n",
                    seqid, gt_error_get(error));
    }
    GtRange seqrange;
    gt_feature_index_get_range_for_seqid(stream->in_loci, &seqrange, seqid,
                                         error);
    if(gt_error_is_set(error))
    {
      gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_stream_"
                    "parse] error retrieving range for sequence '%s': %s\n",
                    seqid, gt_error_get(error));
    }

    if(seqloci == NULL)
      continue;
    GtUword nloci = gt_array_size(seqloci);
    if(nloci > 1)
      gt_array_sort(seqloci, (GtCompare)agn_genome_node_compare);

    // Handle trivial case
    if(nloci == 0)
    {
      AgnLocus *ilocus = agn_locus_new(seqidstr);
      GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
      gt_genome_node_set_range(ilocus, &seqrange);
      gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
      if(gt_error_is_set(error))
      {
        gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_stream_"
                      "parse] error storing feature node in index: %s\n",
                      gt_error_get(error));
      }
      gt_array_delete(seqloci);
      gt_str_delete(seqidstr);
      continue;
    }
    
    // Handle initial ilocus
    AgnLocus *l1 = *(AgnLocus **)gt_array_get(seqloci, 0);
    GtRange r1 = gt_genome_node_get_range(l1);
    if(r1.start >= seqrange.start + (2*delta))
    {
      if(nloci == 1)
      {
        GtRange newrange = { r1.start - delta, r1.end };
        gt_genome_node_set_range(l1, &newrange);
      }
    
      if(!skipterminal)
      {
        AgnLocus *ilocus = agn_locus_new(seqidstr);
        GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
        GtRange newrange = { seqrange.start, r1.start - delta - 1 };
        gt_genome_node_set_range(ilocus, &newrange);
        gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
        if(gt_error_is_set(error))
        {
          gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_"
                        "stream_parse] error storing feature node in index: "
                        "%s\n", gt_error_get(error));
        }
      }
    }
    else
    {
      GtRange newrange = { seqrange.start, r1.end };
      gt_genome_node_set_range(l1, &newrange);
    }
    
    // Handle internal iloci
    for(i = 0; i < nloci - 1; i++)
    {
      AgnLocus *llocus = *(AgnLocus **)gt_array_get(seqloci, i);
      AgnLocus *rlocus = *(AgnLocus **)gt_array_get(seqloci, i+1);
      GtRange lrange = gt_genome_node_get_range(llocus);
      GtRange rrange = gt_genome_node_get_range(rlocus);
    
      if(lrange.end + delta >= rrange.start)
      {
        GtRange newlrange = { lrange.start, rrange.start - 1 };
        gt_genome_node_set_range(llocus, &newlrange);
        GtRange newrrange = { lrange.end + 1, rrange.end };
        gt_genome_node_set_range(rlocus, &newrrange);
      }
      else if(lrange.end + (2*delta) >= rrange.start)
      {
        GtRange newlrange = { lrange.start, lrange.end + delta };
        gt_genome_node_set_range(llocus, &newlrange);
        GtRange newrrange = { rrange.start - delta, rrange.end };
        gt_genome_node_set_range(rlocus, &newrrange);
      }
      else if(lrange.end + (3*delta) >= rrange.start)
      {
        GtUword midpoint = (lrange.end + rrange.start) / 2;
        GtRange newlrange = { lrange.start, midpoint };
        gt_genome_node_set_range(llocus, &newlrange);
        GtRange newrrange = { midpoint + 1, rrange.end };
        gt_genome_node_set_range(rlocus, &newrrange);
      }
      else
      {
        GtRange newlrange = { lrange.start, lrange.end + delta };
        gt_genome_node_set_range(llocus, &newlrange);
        GtRange newrrange = { rrange.start - delta, rrange.end };
        gt_genome_node_set_range(rlocus, &newrrange);
        AgnLocus *ilocus = agn_locus_new(seqidstr);
        GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
        GtRange newrange = { lrange.end + delta + 1, rrange.start - delta - 1 };
        gt_genome_node_set_range(ilocus, &newrange);
        gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
        if(gt_error_is_set(error))
        {
          gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_"
                        "stream_parse] error storing feature node in index: "
                        "%s\n", gt_error_get(error));
        }
      }
    } // end handle internal loci
    
    // Handle terminal ilocus
    l1 = *(AgnLocus **)gt_array_get(seqloci, nloci - 1);
    r1 = gt_genome_node_get_range(l1);
    if(r1.end <= seqrange.end - (2*delta))
    {
      GtRange newrange = { r1.start, r1.end + delta };
      gt_genome_node_set_range(l1, &newrange);
    
      if(!skipterminal)
      {
        AgnLocus *ilocus = agn_locus_new(seqidstr);
        GtFeatureNode *ilocusfn = gt_feature_node_cast(ilocus);
        GtRange nrange = { r1.end + delta + 1, seqrange.end };
        gt_genome_node_set_range(ilocus, &nrange);
        gt_feature_index_add_feature_node(stream->out_loci, ilocusfn, error);
        if(gt_error_is_set(error))
        {
          gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_"
                        "stream_parse] error storing feature node in index: "
                        "%s\n", gt_error_get(error));
        }
      }
    }
    else
    {
      GtRange newrange = { r1.start, seqrange.end };
      gt_genome_node_set_range(l1, &newrange);
    }
    if(nloci == 1)
    {
      GtFeatureNode *l1fn = gt_feature_node_cast(l1);
      gt_feature_index_add_feature_node(stream->out_loci, l1fn, error);
      if(gt_error_is_set(error))
      {
        gt_logger_log(stream->logger, "[AgnIntervalLocusStream::ilocus_"
                      "stream_parse] error storing feature node in index: "
                      "%s\n", gt_error_get(error));
      }
    }

    gt_str_delete(seqidstr);
  }

  gt_error_delete(error);
  gt_str_array_delete(seqids);
}
