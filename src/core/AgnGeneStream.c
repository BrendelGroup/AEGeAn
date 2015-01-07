/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
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
  GtNodeStream *ns, *current_stream, *last_stream;
  AgnGeneStream *stream;
  agn_assert(in_stream && logger);

  ns = gt_node_stream_create(gene_stream_class(), false);
  stream = gene_stream_cast(ns);
  stream->logger = logger;
  stream->streams = gt_queue_new();
  gt_queue_add(stream->streams, gt_node_stream_ref(in_stream));
  last_stream = in_stream;

  current_stream = agn_infer_cds_stream_new(last_stream, logger);
  gt_queue_add(stream->streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_infer_exons_stream_new(last_stream, logger);
  gt_queue_add(stream->streams, current_stream);
  last_stream = current_stream;

  GtHashmap *typestokeep = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                          gt_free_func);
  gt_hashmap_add(typestokeep, gt_cstr_dup("gene"), gt_cstr_dup("gene"));
  current_stream = agn_filter_stream_new(last_stream, typestokeep);
  gt_queue_add(stream->streams, current_stream);
  last_stream = current_stream;
  gt_hashmap_delete(typestokeep);

  stream->in_stream = last_stream;
  return ns;
}

bool agn_gene_stream_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  gene_stream_test_data(queue);
  agn_assert(gt_queue_size(queue) == 3);

  GtFeatureNode *fn = gt_queue_get(queue);
  GtArray *mrnas = agn_typecheck_select(fn, agn_typecheck_mrna);
  bool test1 = (gt_array_size(mrnas) == 2);
  if(test1)
  {
    GtGenomeNode **fn1 = gt_array_get(mrnas, 0);
    GtGenomeNode **fn2 = gt_array_get(mrnas, 1);
    GtRange range1 = gt_genome_node_get_range(*fn1);
    GtRange range2 = gt_genome_node_get_range(*fn2);
    test1 = (range1.start == 5000 && range1.end == 9000 &&
             range2.start == 6000 && range2.end == 9000);
  }
  agn_unit_test_result(test, "1 gene, 2 mRNAs", test1);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(mrnas);

  fn = gt_queue_get(queue);
  mrnas = agn_typecheck_select(fn, agn_typecheck_mrna);
  bool test2 = (gt_array_size(mrnas) == 1);
  if(test2)
  {
    GtGenomeNode **fn1 = gt_array_get(mrnas, 0);
    GtRange range1 = gt_genome_node_get_range(*fn1);
    test2 = (range1.start == 15000 && range1.end == 19000);
  }
  agn_unit_test_result(test, "1 gene, 2 mRNAs (one sans CDS)", test2);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(mrnas);

  fn = gt_queue_get(queue);
  mrnas = agn_typecheck_select(fn, agn_typecheck_mrna);
  bool test3 = gt_array_size(mrnas) == 1;
  if(test3)
  {
    GtArray *utr3p = agn_typecheck_select(fn, agn_typecheck_utr3p);
    GtArray *utr5p = agn_typecheck_select(fn, agn_typecheck_utr5p);
    GtRange range = gt_genome_node_get_range((GtGenomeNode *)fn);
    test3 = (gt_array_size(utr3p) == 1 && gt_array_size(utr5p) == 0);
    if(test3)
    {
      GtGenomeNode *utr = *(GtGenomeNode **)gt_array_get(utr3p, 0);
      GtRange utrrange = gt_genome_node_get_range(utr);
      test3 = (utrrange.start == 32926 && utrrange.end == 33035 &&
               range.start == 30600 && range.end == 33035);
    }
    gt_array_delete(utr3p);
    gt_array_delete(utr5p);
  }
  agn_unit_test_result(test, "mRNA only", test3);
  gt_genome_node_delete((GtGenomeNode *)fn);
  gt_array_delete(mrnas);

  gt_queue_delete(queue);
  return agn_unit_test_success(test);
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

    agn_assert(agn_typecheck_gene(fn));

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
  GtError *error = gt_error_new();
  const char *file = "data/gff3/gene-stream-data.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &file);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  FILE *log = fopen("/dev/null", "w");
  if(log == NULL)
  {
    fprintf(stderr, "[AgnGeneStream::gene_stream_test_data] error opening "
            "/dev/null");
  }
  GtLogger *logger = gt_logger_new(true, "", log);
  GtNodeStream *stream = agn_gene_stream_new(gff3in, logger);
  GtArray *feats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(stream, feats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnGeneStream::gene_stream_test_data] error processing "
            "features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(stream);
  gt_node_stream_delete(arraystream);
  gt_logger_delete(logger);
  if(log != NULL)
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
