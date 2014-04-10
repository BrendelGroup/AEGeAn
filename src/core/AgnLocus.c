#include <math.h>
#include <string.h>
#include "core/array_api.h"
#include "extended/feature_node_iterator_api.h"
#include "AgnLocus.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function The Bron-Kerbosch algorithm is an algorithm for enumerating all
 * maximal cliques in an undirected graph. See the `algorithm's Wikipedia entry
 * <http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm>`_
 * for a description of ``R``, ``P``, and ``X``. All maximal cliques will be
 * stored in ``cliques``. If ``skipsimplecliques`` is true, cliques containing a
 * single item will not be stored.
 */
static void locus_bron_kerbosch(GtArray *R, GtArray *P, GtArray *X,
                                GtArray *cliques, AgnSequenceRegion *region,
                                bool skipsimplecliques);

/**
 * @function ``GtFree`` function: treats each entry in the array as an
 * ``AgnCliquePair **``, dereferences & deletes each entry, and deletes the
 * array.
 */
static void locus_clique_array_delete(GtArray *array);

/**
 * @function If reference transcripts belonging to the same locus overlap, they
 * must be separated before comparison with prediction transcript models (and
 * vice versa). This is an instance of the maximal clique enumeration problem
 * (NP-complete), for which the Bron-Kerbosch algorithm provides a solution.
 */
static GtArray *locus_enumerate_cliques(AgnLocus *locus, GtArray *trans);

/**
 * @function Once all reference transcript cliques and prediction transcript
 * cliques have been enumerated, this function enumerates every possible
 * pairing of 1 reference clique and 1 prediction clique.
 */
static GtArray *locus_enumerate_pairs(AgnLocus *locus, GtArray *refrcliques,
                                      GtArray *predcliques);

/**
 * @function Determine which clique pairs will actually be reported.
 */
static void locus_select_pairs(AgnLocus *locus, GtArray *refrcliques,
                               GtArray *predcliques, GtArray *clique_pairs);

/**
 * @function Generate data for unit testing.
 */
static void locus_test_data(GtQueue *queue);

/**
 * @function For a set of features, we can construct a graph where each node
 * represents a feature and where two nodes are connected if the corresponding
 * features do not overlap. This function returns the intersection of
 * ``trans`` with the neighbors of ``gn`` (where a "neighbor" refers to an
 * adjacent node).
 */
static GtArray *locus_transcript_neighbors(GtGenomeNode *gn, GtArray *trans);

/**
 * @function Test whether a transcript should be filtered.
 */
static bool locus_gene_source_test(AgnLocus *locus, GtFeatureNode *transcript,
                                   AgnComparisonSource source);

/**
 * @function Update the locus feature when a transcript is added.
 */
static void locus_update_range(AgnLocus *locus, GtFeatureNode *transcript);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_locus_add(AgnLocus *locus, GtFeatureNode *feature,
                   AgnComparisonSource source)
{
  gt_genome_node_ref((GtGenomeNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)locus, feature);
  locus_update_range(locus, feature);

  if(source == DEFAULTSOURCE)
    return;

  GtHashmap *feats = gt_genome_node_get_user_data(locus, "refrfeats");
  if(feats == NULL)
  {
    feats = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    gt_genome_node_add_user_data(locus, "refrfeats", feats,
                                 (GtFree)gt_hashmap_delete);
  }
  if(source == REFERENCESOURCE)
    gt_hashmap_add(feats, feature, feature);

  feats = gt_genome_node_get_user_data(locus, "predfeats");
  if(feats == NULL)
  {
    feats = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    gt_genome_node_add_user_data(locus, "predfeats", feats,
                                 (GtFree)gt_hashmap_delete);
  }
  if(source == PREDICTIONSOURCE)
    gt_hashmap_add(feats, feature, feature);
}

AgnLocus *agn_locus_clone(AgnLocus *locus)
{
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  AgnLocus *newlocus = agn_locus_new(seqid);
  GtFeatureNode *locusfn = gt_feature_node_cast(locus);
  GtFeatureNode *newlocusfn = gt_feature_node_cast(newlocus);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(locusfn);
  GtFeatureNode *fn;
  for(fn  = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn  = gt_feature_node_iterator_next(iter))
  {
    gt_genome_node_ref((GtGenomeNode *)fn);
    gt_feature_node_add_child(newlocusfn, fn);
    locus_update_range(newlocus, fn);
  }

  GtHashmap *refr_feats = gt_genome_node_get_user_data(locus, "refrfeats");
  if(refr_feats != NULL)
  {
    gt_hashmap_ref(refr_feats);
    gt_genome_node_add_user_data(newlocus, "refrfeats", refr_feats,
                                 (GtFree)gt_hashmap_delete);
  }
  GtHashmap *pred_feats = gt_genome_node_get_user_data(locus, "predfeats");
  if(pred_feats != NULL)
  {
    gt_hashmap_ref(pred_feats);
    gt_genome_node_add_user_data(newlocus, "predfeats", pred_feats,
                                 (GtFree)gt_hashmap_delete);
  }

  GtArray *pairs2report = gt_genome_node_get_user_data(locus, "pairs2report");
  if(pairs2report != NULL)
  {
    gt_array_ref(pairs2report);
    gt_genome_node_add_user_data(newlocus, "pairs2report", pairs2report,
                                 (GtFree)gt_array_delete);
  }

  AgnComparison *compstats = gt_genome_node_get_user_data(locus, "compstats");
  gt_assert(compstats != NULL);
  AgnComparison *newstats = gt_malloc( sizeof(AgnComparison) );
  agn_comparison_init(newstats);
  agn_comparison_aggregate(newstats, compstats);
  gt_genome_node_add_user_data(newlocus, "compstats", compstats,
                               (GtFree)gt_free_func);

  GtArray *uniqrefr = gt_genome_node_get_user_data(locus, "uniqrefr");
  if(uniqrefr != NULL)
  {
    gt_array_ref(uniqrefr);
    gt_genome_node_add_user_data(newlocus, "uniqrefr", uniqrefr,
                                 (GtFree)gt_array_delete);
  }
  GtArray *uniqpred = gt_genome_node_get_user_data(locus, "uniqpred");
  if(uniqpred != NULL)
  {
    gt_array_ref(uniqpred);
    gt_genome_node_add_user_data(newlocus, "uniqpred", uniqpred,
                                 (GtFree)gt_array_delete);
  }

  return newlocus;
}

GtUword agn_locus_cds_length(AgnLocus *locus, AgnComparisonSource src)
{
  GtUword length = 0;
  GtArray *mrnas = agn_locus_mrnas(locus, src);
  while(gt_array_size(mrnas) > 0)
  {
    GtFeatureNode **mrna = gt_array_pop(mrnas);
    length += agn_mrna_cds_length(*mrna);
  }
  gt_array_delete(mrnas);

  return length;
}

void agn_locus_comparative_analysis(AgnLocus *locus, GtUword maxtranscripts,
                                    GtUword maxpairs, GtLogger *logger)
{
  GtArray *pairs2report = gt_genome_node_get_user_data(locus, "pairs2report");
  if(pairs2report != NULL)
    return;

  GtUword numrefr = agn_locus_num_refr_mrnas(locus);
  GtUword numpred = agn_locus_num_pred_mrnas(locus);
  bool noneedanaly = numrefr > maxtranscripts || numpred > maxtranscripts;
  if(maxtranscripts > 0 && noneedanaly)
  {
    GtStr *seqid = gt_genome_node_get_seqid(locus);
    GtRange range = gt_genome_node_get_range(locus);
    gt_logger_log(logger, "warning: locus %s[%lu, %lu] includes %lu reference "
                  "transcripts and %lu prediction transcripts (maximum is %lu);"
                  " skipping this locus\n", gt_str_get(seqid), range.start,
                  range.end, numrefr, numpred, maxtranscripts);
    return;
  }
  GtArray *refr_trans = agn_locus_refr_mrnas(locus);
  GtArray *refrcliques = locus_enumerate_cliques(locus, refr_trans);
  gt_array_delete(refr_trans);
  GtArray *pred_trans = agn_locus_pred_mrnas(locus);
  GtArray *predcliques = locus_enumerate_cliques(locus, pred_trans);
  gt_array_delete(pred_trans);

  GtUword numpairs = gt_array_size(refrcliques) * gt_array_size(predcliques);
  if(maxpairs > 0 && numpairs > maxpairs)
  {
    GtStr *seqid = gt_genome_node_get_seqid(locus);
    GtRange range = gt_genome_node_get_range(locus);
    gt_logger_log(logger, "warning: locus %s[%lu, %lu] includes %lu possible "
                  "transcript clique pairs (maximum is %lu); skipping this "
                  "locus", gt_str_get(seqid), range.start, range.end, numpairs,
                  maxpairs);
    gt_array_delete(refrcliques);
    gt_array_delete(predcliques);
    return;
  }
  GtArray *clique_pairs = locus_enumerate_pairs(locus,refrcliques,predcliques);
  gt_assert(numpairs == gt_array_size(clique_pairs));
  gt_array_sort(clique_pairs, (GtCompare)agn_clique_pair_compare_reverse);
  locus_select_pairs(locus, refrcliques, predcliques, clique_pairs);

  gt_array_delete(refrcliques);
  gt_array_delete(predcliques);
  gt_array_delete(clique_pairs);
}

int agn_locus_array_compare(const void *p1, const void *p2)
{
  AgnLocus *l1 = *(AgnLocus **)p1;
  AgnLocus *l2 = *(AgnLocus **)p2;

  GtStr *seq1 = gt_genome_node_get_seqid(l1);
  GtStr *seq2 = gt_genome_node_get_seqid(l1);
  if(gt_str_cmp(seq1, seq2) != 0)
    return gt_str_cmp(seq1, seq2);

  GtRange l1r = gt_genome_node_get_range(l1);
  GtRange l2r = gt_genome_node_get_range(l2);
  return gt_range_compare(&l1r, &l2r);
}

void agn_locus_comparison_aggregate(AgnLocus *locus, AgnComparison *comp)
{
  AgnComparison *stats = gt_genome_node_get_user_data(locus, "compstats");
  gt_assert(stats != NULL);
  agn_comparison_aggregate(comp, stats);
}

void agn_locus_delete(AgnLocus *locus)
{
  gt_genome_node_delete(locus);
}

GtUword agn_locus_exon_num(AgnLocus *locus, AgnComparisonSource src)
{
  GtUword count = 0;
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_gene(feature));
    if(locus_gene_source_test(locus, feature, src) == false)
      continue;

    GtFeatureNodeIterator *subiter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *subfeature;
    for(subfeature  = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature  = gt_feature_node_iterator_next(subiter))
    {
      if(agn_typecheck_exon(subfeature))
        count++;
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  return count;
}

bool agn_locus_filter_test(AgnLocus *locus, AgnLocusFilter *filter,
                           AgnComparisonSource src)
{
  GtUword value = filter->function(locus, src);
  switch(filter->operator)
  {
    case AGN_LOCUS_FILTER_EQ:
      if(value == filter->testvalue) return true;
      break;
    case AGN_LOCUS_FILTER_NE:
      if(value != filter->testvalue) return true;
      break;
    case AGN_LOCUS_FILTER_GT:
      if(value > filter->testvalue) return true;
      break;
    case AGN_LOCUS_FILTER_GE:
      if(value >= filter->testvalue) return true;
      break;
    case AGN_LOCUS_FILTER_LT:
      if(value < filter->testvalue) return true;
      break;
    case AGN_LOCUS_FILTER_LE:
      if(value <= filter->testvalue) return true;
      break;
  }
  return false;
}

GtArray *agn_locus_get_unique_pred_cliques(AgnLocus *locus)
{
  return gt_genome_node_get_user_data(locus, "unique_pred");
}

GtArray *agn_locus_get_unique_refr_cliques(AgnLocus *locus)
{
  return gt_genome_node_get_user_data(locus, "unique_refr");
}

GtArray *agn_locus_gene_ids(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *ids = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    bool isgene = agn_typecheck_gene(feature);
    bool meets_crit = locus_gene_source_test(locus, feature, src);
    if(isgene && meets_crit)
    {
      const char *id = gt_feature_node_get_attribute(feature, "ID");
      gt_array_add(ids, id);
    }
  }
  gt_feature_node_iterator_delete(iter);

  return ids;
}

GtArray *agn_locus_mrnas(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *mrnas = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    bool isgene = agn_typecheck_gene(feature);
    bool meets_crit = locus_gene_source_test(locus, feature, src);
    if(!isgene || !meets_crit)
      continue;

    GtFeatureNodeIterator *subiter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *subfeature;
    for(subfeature  = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature  = gt_feature_node_iterator_next(subiter))
    {
      if(agn_typecheck_mrna(subfeature))
        gt_array_add(mrnas, subfeature);
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  return mrnas;
}

GtArray *agn_locus_mrna_ids(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *ids = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    bool isgene = agn_typecheck_gene(feature);
    bool meets_crit = locus_gene_source_test(locus, feature, src);
    if(!isgene || !meets_crit)
      continue;

    GtFeatureNodeIterator *subiter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *subfeature;
    for(subfeature  = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature  = gt_feature_node_iterator_next(subiter))
    {
      if(agn_typecheck_mrna(subfeature))
      {
        const char *id = gt_feature_node_get_attribute(feature, "ID");
        gt_array_add(ids, id);
      }
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  return ids;
}

GtUword agn_locus_mrna_num(AgnLocus *locus, AgnComparisonSource src)
{
  GtUword count = 0;
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    bool isgene = agn_typecheck_gene(feature);
    bool meets_crit = locus_gene_source_test(locus, feature, src);
    if(!isgene || !meets_crit)
      continue;

    GtFeatureNodeIterator *subiter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *subfeature;
    for(subfeature  = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature  = gt_feature_node_iterator_next(subiter))
    {
      if(agn_typecheck_mrna(subfeature))
        count++;
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  return count;
}

AgnLocus *agn_locus_new(GtStr *seqid)
{
  AgnLocus *locus = gt_feature_node_new(seqid, "locus", 0, 0, GT_STRAND_BOTH);
  AgnComparison *compstats = gt_malloc( sizeof(AgnComparison) );
  agn_comparison_init(compstats);
  gt_genome_node_add_user_data(locus, "compstats", compstats,
                               (GtFree)gt_free_func);
  return locus;
}

GtArray *agn_locus_pairs_to_report(AgnLocus *locus)
{
  return gt_genome_node_get_user_data(locus, "pairs2report");
}

#ifndef WITHOUT_CAIRO
void agn_locus_png_track_selector(GtBlock *block, GtStr *track, void *data)
{
  GtFeatureNode *fn = gt_block_get_top_level_feature(block);
  GtGenomeNode *gn = (GtGenomeNode *)fn;
  const char *filename = gt_genome_node_get_filename(gn);
  AgnLocusPngMetadata *metadata = data;
  char trackname[512];

  if(strcmp(filename, metadata->refrfile) == 0)
  {
    if(strcmp(metadata->refrlabel, "") == 0)
    {
      sprintf(trackname, "Reference annotations (%s)", metadata->refrfile);
      gt_str_set(track, trackname);
    }
    else
    {
      sprintf(trackname, "%s (Reference)", metadata->refrlabel);
      gt_str_set(track, trackname);
    }
    return;
  }
  else if(strcmp(filename, metadata->predfile) == 0)
  {
    if(strcmp(metadata->predlabel, "") == 0)
    {
      sprintf(trackname, "Prediction annotations (%s)", metadata->predfile);
      gt_str_set(track, trackname);
    }
    else
    {
      sprintf(trackname, "%s (Prediction)", metadata->predlabel);
      gt_str_set(track, trackname);
    }
    return;
  }

  fprintf(stderr, "Error: unknown filename '%s'\n", filename);
  exit(1);
}
#endif

#ifndef WITHOUT_CAIRO
void agn_locus_print_png(AgnLocus *locus, AgnLocusPngMetadata *metadata)
{
  GtError *error = gt_error_new();
  GtFeatureIndex *index = gt_feature_index_memory_new();
  GtUword i;

  GtArray *refr_trans = agn_locus_refr_mrnas(locus);
  for(i = 0; i < gt_array_size(refr_trans); i++)
  {
    GtFeatureNode *trans = *(GtFeatureNode **)gt_array_get(refr_trans, i);
    gt_feature_index_add_feature_node(index, trans, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "error: %s\n", gt_error_get(error));
      exit(1);
    }
  }
  gt_array_delete(refr_trans);

  GtArray *pred_trans = agn_locus_pred_mrnas(locus);
  for(i = 0; i < gt_array_size(pred_trans); i++)
  {
    GtFeatureNode *trans = *(GtFeatureNode **)gt_array_get(pred_trans, i);
    gt_feature_index_add_feature_node(index, trans, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "error: %s\n", gt_error_get(error));
      exit(1);
    }
  }
  gt_array_delete(pred_trans);

  // Generate the graphic...this is going to get a bit hairy
  GtStyle *style;
  if(!(style = gt_style_new(error)))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_style_load_file(style, metadata->stylefile, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange locusrange = gt_genome_node_get_range(locus);
  GtDiagram *diagram = gt_diagram_new(index, gt_str_get(seqid), &locusrange,
                                      style, error);
  gt_diagram_set_track_selector_func(
                              diagram,
                              (GtTrackSelectorFunc)agn_locus_png_track_selector,
                              metadata
  );
  GtLayout *layout = gt_layout_new(diagram, metadata->graphic_width, style,
                                   error);
  if(!layout)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  gt_layout_set_track_ordering_func(
                                layout,
                                (GtTrackOrderingFunc)metadata->track_order_func,
                                NULL
  );
  GtUword image_height;
  if(gt_layout_get_height(layout, &image_height, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  GtCanvas *canvas = gt_canvas_cairo_file_new(style, GT_GRAPHICS_PNG,
                                              metadata->graphic_width,
                                              image_height, NULL, error);
  if(!canvas)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_layout_sketch(layout, canvas, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }
  if(gt_canvas_cairo_file_to_file((GtCanvasCairoFile*) canvas,
                                  metadata->filename, error))
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
    exit(EXIT_FAILURE);
  }

  gt_feature_index_delete(index);
  gt_canvas_delete(canvas);
  gt_layout_delete(layout);
  gt_diagram_delete(diagram);
  gt_style_delete(style);
  gt_error_delete(error);
}
#endif

void agn_locus_print_transcript_mapping(AgnLocus *locus, FILE *outstream)
{
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange locusrange = gt_genome_node_get_range(locus);
  GtArray *transids = agn_locus_get_mrna_ids(locus);
  while(gt_array_size(transids) > 0)
  {
    const char **transid = gt_array_pop(transids);
    fprintf(outstream, "%s\t%s:%lu-%lu\n", *transid, gt_str_get(seqid),
            locusrange.start, locusrange.end);
  }
  gt_array_delete(transids);
}

void agn_locus_set_range(AgnLocus *locus, GtUword start, GtUword end)
{
  GtRange range = { start, end };
  gt_genome_node_set_range(locus, &range);
}

double agn_locus_splice_complexity(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *trans = agn_locus_mrnas(locus, src);
  double sc = agn_calc_splice_complexity(trans);
  gt_array_delete(trans);
  return sc;
}

bool agn_locus_unit_test(AgnUnitTest *test)
{
  GtQueue *queue = gt_queue_new();
  locus_test_data(queue);
  gt_assert(gt_queue_size(queue) == 4);

  GtLogger *logger = gt_logger_new(true, "", stderr);

  AgnComparison c;
  c.cds_nuc_stats.tp = 870;
  c.cds_nuc_stats.fn = 0;
  c.cds_nuc_stats.fp = 0;
  c.cds_nuc_stats.tn = 2240;
  c.cds_nuc_stats.mc = 1.0;
  c.cds_nuc_stats.cc = 1.0;
  c.cds_nuc_stats.sn = 1.0;
  c.cds_nuc_stats.sp = 1.0;
  c.cds_nuc_stats.f1 = 1.0;
  c.cds_nuc_stats.ed = 0.0;
  c.utr_nuc_stats.tp = 337;
  c.utr_nuc_stats.fn = 11;
  c.utr_nuc_stats.fp = 131;
  c.utr_nuc_stats.tn = 2631;
  c.utr_nuc_stats.mc = 0.954341;
  c.utr_nuc_stats.cc = 0.811995;
  c.utr_nuc_stats.sn = 0.968391;
  c.utr_nuc_stats.sp = 0.720085;
  c.utr_nuc_stats.f1 = 0.825980;
  c.utr_nuc_stats.ed = 0.155762;
  c.cds_struc_stats.correct = 12;
  c.cds_struc_stats.missing = 0;
  c.cds_struc_stats.wrong   = 0;
  c.cds_struc_stats.sn = 1.0;
  c.cds_struc_stats.sp = 1.0;
  c.cds_struc_stats.f1 = 1.0;
  c.cds_struc_stats.ed = 0.0;
  c.exon_struc_stats.correct = 10;
  c.exon_struc_stats.missing = 2;
  c.exon_struc_stats.wrong   = 2;
  c.exon_struc_stats.sn = 0.833333;
  c.exon_struc_stats.sp = 0.833333;
  c.exon_struc_stats.f1 = 0.833333;
  c.exon_struc_stats.ed = 0.166667;
  c.utr_struc_stats.correct = 0;
  c.utr_struc_stats.missing = 1;
  c.utr_struc_stats.wrong   = 2;
  c.utr_struc_stats.sn = 0.0;
  c.utr_struc_stats.sp = 0.0;
  c.utr_struc_stats.f1 = 0 / (c.utr_nuc_stats.f1-c.utr_nuc_stats.f1);
  c.utr_struc_stats.ed = 1.0;
  c.overall_matches = 2968;
  c.overall_length  = 3110;

  AgnLocus *locus = gt_queue_get(queue);
  agn_locus_comparative_analysis(locus, 0, 0, logger);
  AgnComparison stats;
  agn_comparison_init(&stats);
  agn_locus_comparison_aggregate(locus, &stats);
  agn_comparison_resolve(&stats);
  bool grapetest1 = agn_comparison_test(&stats, &c);
  agn_unit_test_result(test, "grape test 1", grapetest1);
  agn_locus_delete(locus);


  c.cds_nuc_stats.tp = 379;
  c.cds_nuc_stats.fn = 194;
  c.cds_nuc_stats.fp = 506;
  c.cds_nuc_stats.tn = 3147;
  c.cds_nuc_stats.mc = 0.834359;
  c.cds_nuc_stats.cc = 0.439970;
  c.cds_nuc_stats.sn = 0.661431;
  c.cds_nuc_stats.sp = 0.428249;
  c.cds_nuc_stats.f1 = 0.519890;
  c.cds_nuc_stats.ed = 0.455160;
  c.utr_nuc_stats.tp = 102;
  c.utr_nuc_stats.fn = 0;
  c.utr_nuc_stats.fp = 114;
  c.utr_nuc_stats.tn = 4010;
  c.utr_nuc_stats.mc = 0.973024;
  c.utr_nuc_stats.cc = 0.677620;
  c.utr_nuc_stats.sn = 1.0;
  c.utr_nuc_stats.sp = 0.472222;
  c.utr_nuc_stats.f1 = 0.641509;
  c.utr_nuc_stats.ed = 0.263889;
  c.cds_struc_stats.correct = 0;
  c.cds_struc_stats.missing = 4;
  c.cds_struc_stats.wrong   = 1;
  c.cds_struc_stats.sn = 0.0;
  c.cds_struc_stats.sp = 0.0;
  c.cds_struc_stats.f1 = 0 / (c.utr_nuc_stats.f1-c.utr_nuc_stats.f1);
  c.cds_struc_stats.ed = 1.0;
  c.exon_struc_stats.correct = 2;
  c.exon_struc_stats.missing = 3;
  c.exon_struc_stats.wrong   = 1;
  c.exon_struc_stats.sn = 0.400000;
  c.exon_struc_stats.sp = 0.666667;
  c.exon_struc_stats.f1 = 0.500000;
  c.exon_struc_stats.ed = 0.466667;
  c.utr_struc_stats.correct = 1;
  c.utr_struc_stats.missing = 1;
  c.utr_struc_stats.wrong   = 2;
  c.utr_struc_stats.sn = 0.500000;
  c.utr_struc_stats.sp = 0.333333;
  c.utr_struc_stats.f1 = 0.400000;
  c.utr_struc_stats.ed = 0.583333;
  c.overall_matches = 3379;
  c.overall_length  = 4226;

  locus = gt_queue_get(queue);
  agn_locus_comparative_analysis(locus, 0, 0, logger);
  agn_comparison_init(&stats);
  agn_locus_comparison_aggregate(locus, &stats);
  agn_comparison_resolve(&stats);
  bool grapetest2 = agn_comparison_test(&stats, &c);
  agn_unit_test_result(test, "grape test 2", grapetest2);
  agn_locus_delete(locus);

  AgnLocus *locus1 = gt_queue_get(queue);
  AgnLocus *locus2 = gt_queue_get(queue);
  AgnLocusFilter filter = { agn_locus_cds_length, 699, AGN_LOCUS_FILTER_LT };
  bool cdslengthtest = !agn_locus_filter_test(locus1, &filter, DEFAULTSOURCE) &&
                       !agn_locus_filter_test(locus2, &filter, DEFAULTSOURCE);
  filter.testvalue = 750;
  cdslengthtest = cdslengthtest &&
                  agn_locus_filter_test(locus1, &filter, DEFAULTSOURCE) &&
                  !agn_locus_filter_test(locus2, &filter, DEFAULTSOURCE);
  filter.testvalue = 813;
  filter.operator = AGN_LOCUS_FILTER_LE;
  cdslengthtest = cdslengthtest &&
                  agn_locus_filter_test(locus1, &filter, DEFAULTSOURCE) &&
                  agn_locus_filter_test(locus2, &filter, DEFAULTSOURCE);
  agn_unit_test_result(test, "filter by CDS length", cdslengthtest);

  filter.function = agn_locus_exon_num;
  filter.testvalue = 3;
  filter.operator = AGN_LOCUS_FILTER_NE;
  bool exonnumtest = !agn_locus_filter_test(locus1, &filter, DEFAULTSOURCE) &&
                      agn_locus_filter_test(locus2, &filter, DEFAULTSOURCE);
  filter.testvalue = 7;
  exonnumtest = exonnumtest &&
                agn_locus_filter_test(locus1, &filter, DEFAULTSOURCE) &&
                !agn_locus_filter_test(locus2, &filter, DEFAULTSOURCE);
  agn_unit_test_result(test, "filter by exon number", exonnumtest);
  agn_locus_delete(locus1);
  agn_locus_delete(locus2);

  gt_logger_delete(logger);
  gt_queue_delete(queue);
  return agn_unit_test_success(test);
}

static void locus_bron_kerbosch(GtArray *R, GtArray *P, GtArray *X,
                                GtArray *cliques, AgnSequenceRegion *region,
                                bool skipsimplecliques)
{
  gt_assert(R != NULL && P != NULL && X != NULL && cliques != NULL);

  if(gt_array_size(P) == 0 && gt_array_size(X) == 0)
  {
    if(skipsimplecliques == false || gt_array_size(R) != 1)
    {
      GtUword i;
      AgnTranscriptClique *clique = agn_transcript_clique_new(region);
      for(i = 0; i < gt_array_size(R); i++)
      {
        GtFeatureNode *transcript = *(GtFeatureNode **)gt_array_get(R, i);
        agn_transcript_clique_add(clique, transcript);
      }
      gt_array_add(cliques, clique);
    }
  }

  while(gt_array_size(P) > 0)
  {
    GtGenomeNode *v = *(GtGenomeNode **)gt_array_get(P, 0);

    // newR = R \union {v}
    GtArray *newR = agn_array_copy(R, sizeof(GtGenomeNode *));
    gt_array_add(newR, v);
    // newP = P \intersect N(v)
    GtArray *newP = locus_transcript_neighbors(v, P);
    // newX = X \intersect N(v)
    GtArray *newX = locus_transcript_neighbors(v, X);

    // Recursive call
    // locus_bron_kerbosch(R \union {v}, P \intersect N(v), X \intersect N(X))
    locus_bron_kerbosch(newR, newP, newX, cliques, region, skipsimplecliques);

    // Delete temporary arrays just created
    gt_array_delete(newR);
    gt_array_delete(newP);
    gt_array_delete(newX);

    // P := P \ {v}
    gt_array_rem(P, 0);

    // X := X \union {v}
    gt_array_add(X, v);
  }
}

static void locus_clique_array_delete(GtArray *array)
{
  gt_assert(array != NULL);
  while(gt_array_size(array) > 0)
  {
    AgnCliquePair **pair = gt_array_pop(array);
    agn_clique_pair_delete(*pair);
  }
  gt_array_delete(array);
}

static GtArray *locus_enumerate_cliques(AgnLocus *locus, GtArray *trans)
{
  GtArray *cliques = gt_array_new( sizeof(AgnTranscriptClique *) );
  GtUword numtrans = gt_array_size(trans);
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  GtRange range = gt_genome_node_get_range(locus);
  AgnSequenceRegion region = { seqid, range };

  if(numtrans == 1)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(trans, 0);
    AgnTranscriptClique *clique = agn_transcript_clique_new(&region);
    agn_transcript_clique_add(clique, fn);
    gt_array_add(cliques, clique);
  }
  else
  {
    // First add each transcript as a clique, even if it is not a maximal clique
    GtUword i;
    for(i = 0; i < numtrans; i++)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(trans, i);
      AgnTranscriptClique *clique = agn_transcript_clique_new(&region);
      agn_transcript_clique_add(clique, fn);
      gt_array_add(cliques, clique);
    }

    // Then use the Bron-Kerbosch algorithm to find all maximal cliques
    // containing >1 transcript
    GtArray *R = gt_array_new( sizeof(GtGenomeNode *) );
    GtArray *P = agn_array_copy(trans, sizeof(GtGenomeNode *));
    GtArray *X = gt_array_new( sizeof(GtGenomeNode *) );

    // Initial call: locus_bron_kerbosch(\emptyset, vertex_set, \emptyset )
    locus_bron_kerbosch(R, P, X, cliques, &region, true);

    gt_array_delete(R);
    gt_array_delete(P);
    gt_array_delete(X);
  }

  return cliques;
}

static GtArray *locus_enumerate_pairs(AgnLocus *locus, GtArray *refrcliques,
                                      GtArray *predcliques)
{
  gt_assert(refrcliques != NULL && predcliques != NULL);

  GtArray *clique_pairs = gt_array_new( sizeof(AgnCliquePair *) );
  GtUword i,j;
  for(i = 0; i < gt_array_size(refrcliques); i++)
  {
    AgnTranscriptClique *refr_clique, *pred_clique;
    refr_clique = *(AgnTranscriptClique **)gt_array_get(refrcliques, i);
    for(j = 0; j < gt_array_size(predcliques); j++)
    {
      pred_clique = *(AgnTranscriptClique**)gt_array_get(predcliques, j);
      AgnCliquePair *pair = agn_clique_pair_new(refr_clique, pred_clique);
      gt_array_add(clique_pairs, pair);
    }
  }

  return clique_pairs;
}

static void locus_select_pairs(AgnLocus *locus, GtArray *refrcliques,
                               GtArray *predcliques, GtArray *clique_pairs)
{
  GtHashmap *refrcliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  GtHashmap *predcliques_acctd = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);

  AgnComparison *stats = gt_genome_node_get_user_data(locus, "compstats");
  gt_assert(stats != NULL);
  GtArray *pairs2report = gt_array_new( sizeof(AgnCliquePair *) );
  GtUword i;
  for(i = 0; i < gt_array_size(clique_pairs); i++)
  {
    AgnCliquePair **pair = gt_array_get(clique_pairs, i);
    AgnTranscriptClique *rclique = agn_clique_pair_get_refr_clique(*pair);
    AgnTranscriptClique *pclique = agn_clique_pair_get_pred_clique(*pair);
    if(agn_transcript_clique_has_id_in_hash(rclique, refrcliques_acctd) &&
       agn_transcript_clique_has_id_in_hash(pclique, predcliques_acctd))
    {
      agn_clique_pair_delete(*pair);
    }
    else
    {
      gt_array_add(pairs2report, *pair);
      agn_clique_pair_comparison_aggregate(*pair, stats);
      agn_transcript_clique_put_ids_in_hash(rclique, refrcliques_acctd);
      agn_transcript_clique_put_ids_in_hash(pclique, predcliques_acctd);
    }
  }
  gt_genome_node_add_user_data(locus,"pairs2report",gt_array_ref(pairs2report),
                               (GtFree)locus_clique_array_delete);
  gt_array_delete(pairs2report);
  agn_comparison_resolve(stats);

  GtArray *uniqrefr = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(refrcliques); i++)
  {
    AgnTranscriptClique *refr_clique;
    refr_clique = *(AgnTranscriptClique **)gt_array_get(refrcliques, i);
    if(!agn_transcript_clique_has_id_in_hash(refr_clique, refrcliques_acctd))
    {
      gt_genome_node_ref(refr_clique);
      gt_array_add(uniqrefr, refr_clique);
      agn_transcript_clique_put_ids_in_hash(refr_clique, refrcliques_acctd);
    }
    agn_transcript_clique_delete(refr_clique);
  }
  if(gt_array_size(uniqrefr) > 0)
  {
    gt_genome_node_add_user_data(locus, "uniqrefr", gt_array_ref(uniqrefr),
                                 (GtFree)gt_array_delete);
  }
  gt_array_delete(uniqrefr);

  GtArray *uniqpred = gt_array_new( sizeof(AgnTranscriptClique *) );
  for(i = 0; i < gt_array_size(predcliques); i++)
  {
    AgnTranscriptClique *pred_clique;
    pred_clique = *(AgnTranscriptClique **)gt_array_get(predcliques, i);
    if(!agn_transcript_clique_has_id_in_hash(pred_clique, predcliques_acctd))
    {
      gt_array_add_elem(uniqpred, gt_genome_node_ref(pred_clique),
                        sizeof(GtArray *));
      agn_transcript_clique_put_ids_in_hash(pred_clique, predcliques_acctd);
    }
    agn_transcript_clique_delete(pred_clique);
  }
  if(gt_array_size(uniqpred) > 0)
  {
    gt_genome_node_add_user_data(locus, "uniqpred", gt_array_ref(uniqpred),
                                 (GtFree)gt_array_delete);
  }
  gt_array_delete(uniqpred);

  gt_hashmap_delete(refrcliques_acctd);
  gt_hashmap_delete(predcliques_acctd);
}

static void locus_test_data(GtQueue *queue)
{
  gt_assert(queue != NULL);
  GtArray *refrfeats, *predfeats;

  GtError *error = gt_error_new();
  const char *refrfile = "data/gff3/grape-refr.gff3";
  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  refrfeats = gt_array_new( sizeof(GtFeatureNode *) );
  GtNodeStream *arraystream = gt_array_out_stream_new(gff3in, refrfeats, error);
  int pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnLocus::locus_test_data] error processing reference "
            "features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_array_sort(refrfeats, (GtCompare)agn_genome_node_compare);

  const char *predfile = "data/gff3/grape-pred.gff3";
  gff3in = gt_gff3_in_stream_new_unsorted(1, &predfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  predfeats = gt_array_new( sizeof(GtFeatureNode *) );
  arraystream = gt_array_out_stream_new(gff3in, predfeats, error);
  pullresult = gt_node_stream_pull(arraystream, error);
  if(pullresult == -1)
  {
    fprintf(stderr, "[AgnLocus::locus_test_data] error processing prediction "
            "features: %s\n", gt_error_get(error));
  }
  gt_node_stream_delete(gff3in);
  gt_node_stream_delete(arraystream);
  gt_array_sort(predfeats, (GtCompare)agn_genome_node_compare);

  gt_assert(gt_array_size(refrfeats) == 12 && gt_array_size(predfeats) == 13);

  GtStr *seqid = gt_str_new_cstr("chr8");
  GtFeatureNode *refr = *(GtFeatureNode **)gt_array_get(refrfeats, 2);
  GtFeatureNode *pred = *(GtFeatureNode **)gt_array_get(predfeats, 3);
  AgnLocus *locus = agn_locus_new(seqid);
  agn_locus_add_refr_feature(locus, refr);
  agn_locus_add_pred_feature(locus, pred);
  gt_queue_add(queue, locus);

  refr = *(GtFeatureNode **)gt_array_get(refrfeats, 9);
  pred = *(GtFeatureNode **)gt_array_get(predfeats, 11);
  locus = agn_locus_new(seqid);
  agn_locus_add_refr_feature(locus, refr);
  agn_locus_add_pred_feature(locus, pred);
  gt_queue_add(queue, locus);

  refr = *(GtFeatureNode **)gt_array_get(refrfeats, 0);
  locus = agn_locus_new(seqid);
  agn_locus_add_feature(locus, refr);
  gt_queue_add(queue, locus);

  refr = *(GtFeatureNode **)gt_array_get(refrfeats, 3);
  locus = agn_locus_new(seqid);
  agn_locus_add_feature(locus, refr);
  gt_queue_add(queue, locus);

  while(gt_array_size(refrfeats) > 0)
  {
    GtGenomeNode **gn = gt_array_pop(refrfeats);
    gt_genome_node_delete(*gn);
  }
  gt_array_delete(refrfeats);
  while(gt_array_size(predfeats) > 0)
  {
    GtGenomeNode **gn = gt_array_pop(predfeats);
    gt_genome_node_delete(*gn);
  }
  gt_array_delete(predfeats);

  gt_str_delete(seqid);
  gt_error_delete(error);
}

static GtArray *locus_transcript_neighbors(GtGenomeNode *gn, GtArray *trans)
{
  GtArray *neighbors = gt_array_new( sizeof(GtGenomeNode *) );
  GtUword i;
  for(i = 0; i < gt_array_size(trans); i++)
  {
    GtGenomeNode *other = *(GtGenomeNode **)gt_array_get(trans, i);
    if(other != gn)
    {
      GtRange gn_range = gt_genome_node_get_range(gn);
      GtRange other_range = gt_genome_node_get_range(other);
      if(gt_range_overlap(&gn_range, &other_range) == false)
        gt_array_add(neighbors, other);
    }
  }
  return neighbors;
}

static bool locus_gene_source_test(AgnLocus *locus, GtFeatureNode *transcript,
                                   AgnComparisonSource source)
{
  if(source == DEFAULTSOURCE)
    return true;

  GtHashmap *refr_feats = gt_genome_node_get_user_data(locus, "refrfeats");
  GtHashmap *pred_feats = gt_genome_node_get_user_data(locus, "predfeats");
  bool inrefr = (gt_hashmap_get(refr_feats, transcript) != NULL);
  bool inpred = (gt_hashmap_get(pred_feats, transcript) != NULL);
  if(source == REFERENCESOURCE  && inrefr)
    return true;
  else if(source == PREDICTIONSOURCE && inpred)
    return true;
    
  return false;
}

static void locus_update_range(AgnLocus *locus, GtFeatureNode *transcript)
{
  GtRange locusrange = gt_genome_node_get_range(locus);
  GtRange transrange = gt_genome_node_get_range((GtGenomeNode *)transcript);
  if(locusrange.start == 0 && locusrange.end == 0)
  {
    gt_genome_node_set_range(locus, &transrange);
    return;
  }

  GtRange newrange = gt_range_join(&locusrange, &transrange);
  gt_genome_node_set_range(locus, &newrange);
}
