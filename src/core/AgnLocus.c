#include <string.h>
#include "extended/feature_node_iterator_api.h"
#include "AgnLocus.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Test whether a transcript should be included.
 */
static bool locus_transcript_source_test(AgnLocus *locus,
                                         GtFeatureNode *transcript,
                                         AgnComparisonSource source);

/**
 * @function Update the locus feature when a transcript is added.
 */
static void locus_update_range(AgnLocus *locus, GtFeatureNode *transcript);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_locus_add(AgnLocus *locus, GtFeatureNode *transcript,
                   AgnComparisonSource source)
{
  gt_genome_node_ref((GtGenomeNode *)transcript);
  gt_feature_node_add_child((GtFeatureNode *)locus, transcript);
  locus_update_range(locus, transcript);
  
  const char *key = "refrtrans";
  if(source == DEFAULTSOURCE)
    return;
  else if(source == PREDICTIONSOURCE)
    key = "predtrans";

  GtHashmap *trans = gt_genome_node_get_user_data(locus, key);
  if(trans == NULL)
  {
    trans = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
    gt_genome_node_add_user_data(locus, key, trans, (GtFree)gt_hashmap_delete);
  }
  gt_hashmap_add(trans, transcript, transcript);
}

AgnLocus *agn_locus_clone(AgnLocus *locus)
{
  GtStr *seqid = gt_genome_node_get_seqid(locus);
  AgnLocus *newlocus = agn_locus_new(gt_str_get(seqid));
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

  GtHashmap *refr_trans = gt_genome_node_get_user_data(locus, "refrtrans");
  if(refr_trans != NULL)
  {
    gt_hashmap_ref(refr_trans);
    gt_genome_node_add_user_data(newlocus, "refrtrans", refr_trans,
                                 (GtFree)gt_hashmap_delete);
  }
  GtHashmap *pred_trans = gt_genome_node_get_user_data(locus, "predtrans");
  if(pred_trans != NULL)
  {
    gt_hashmap_ref(pred_trans);
    gt_genome_node_add_user_data(newlocus, "predtrans", pred_trans,
                                 (GtFree)gt_hashmap_delete);
  }

  GtArray *pairs2report = gt_genome_node_get_user_data(locus, "pairs2report");
  if(pairs2report != NULL)
  {
    gt_array_ref(pairs2report);
    gt_genome_node_add_user_data(newlocus, "pairs2report", pairs2report,
                                 (GtFree)gt_array_delete);
  }
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
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(feature));
    if(locus_transcript_source_test(locus, feature, src) == false)
      continue;

    GtFeatureNodeIterator *subiter = gt_feature_node_iterator_new(feature);
    GtFeatureNode *subfeature;
    for(subfeature  = gt_feature_node_iterator_next(subiter);
        subfeature != NULL;
        subfeature  = gt_feature_node_iterator_next(subiter))
    {
      if(agn_typecheck_cds(subfeature))
        length += gt_genome_node_get_length((GtGenomeNode *)subfeature);
    }
    gt_feature_node_iterator_delete(subiter);
  }
  gt_feature_node_iterator_delete(iter);

  return length;
}

void agn_locus_comparative_analysis(AgnLocus *locus)
{
  
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
  if(stats == NULL)
    return;
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
    gt_assert(agn_typecheck_transcript(feature));
    if(locus_transcript_source_test(locus, feature, src) == false)
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

GtArray *agn_locus_get_unique_pred_cliques(AgnLocus *locus)
{
  return gt_genome_node_get_user_data(locus, "unique_pred");
}

GtArray *agn_locus_get_unique_refr_cliques(AgnLocus *locus)
{
  return gt_genome_node_get_user_data(locus, "unique_refr");
}

AgnLocus *agn_locus_new(const char *seqid)
{
  return gt_feature_node_new(gt_str_new_cstr(seqid), "locus", 0, 0,
                             GT_STRAND_BOTH);
}

GtUword agn_locus_num_clique_pairs(AgnLocus *locus)
{
  GtArray *pairs2report = gt_genome_node_get_user_data(locus, "pairs2report");
  return gt_array_size(pairs2report);
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

  GtArray *refr_trans = agn_locus_refr_transcripts(locus);
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

  GtArray *pred_trans = agn_locus_pred_transcripts(locus);
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
  GtArray *transids = agn_locus_get_transcript_ids(locus);
  while(gt_array_size(transids) > 0)
  {
    const char **transid = gt_array_pop(transids);
    fprintf(outstream, "%s\t%s:%lu-%lu\n", *transid, gt_str_get(seqid),
            locusrange.start, locusrange.end);
  }
  gt_array_delete(transids);
}

double agn_locus_splice_complexity(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *trans = agn_locus_transcripts(locus, src);
  double sc = agn_calc_splice_complexity(trans);
  gt_array_delete(trans);
  return sc;
}

GtArray *agn_locus_transcripts(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *transcripts = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(feature));
    if(locus_transcript_source_test(locus, feature, src))
      gt_array_add(transcripts, feature);
  }
  gt_feature_node_iterator_delete(iter);

  return transcripts;
}

GtArray *agn_locus_transcript_ids(AgnLocus *locus, AgnComparisonSource src)
{
  GtArray *ids = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(feature));
    if(locus_transcript_source_test(locus, feature, src))
    {
      const char *id = gt_feature_node_get_attribute(feature, "ID");
      gt_array_add(ids, id);
    }
  }
  gt_feature_node_iterator_delete(iter);

  return ids;
}

GtUword agn_locus_transcript_num(AgnLocus *locus, AgnComparisonSource src)
{
  GtUword count = 0;
  GtFeatureNode *fn = gt_feature_node_cast(locus);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(feature));
    if(locus_transcript_source_test(locus, feature, src))
      count++;
  }
  gt_feature_node_iterator_delete(iter);

  return count;
}

bool agn_locus_unit_test(AgnUnitTest *test)
{
  return false;
}

static bool locus_transcript_source_test(AgnLocus *locus,
                                         GtFeatureNode *transcript,
                                         AgnComparisonSource source)
{
  if(source == DEFAULTSOURCE)
    return true;

  GtHashmap *refr_trans = gt_genome_node_get_user_data(locus, "refrtrans");
  GtHashmap *pred_trans = gt_genome_node_get_user_data(locus, "predtrans");
  bool inrefr = (gt_hashmap_get(refr_trans, transcript) != NULL);
  bool inpred = (gt_hashmap_get(pred_trans, transcript) != NULL);
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
