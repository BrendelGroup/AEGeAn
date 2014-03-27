#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

typedef struct
{
  GtFeatureNode *mrna;
  GtArray *exons;
  GtIntervalTree *exon_tree;
} MrnaStructure;

static MrnaStructure *mrna_structure_new(GtFeatureNode *mrna)
{
  GtUword i;
  MrnaStructure *s = gt_malloc(sizeof(MrnaStructure));
  s->mrna = mrna;
  s->exons = agn_typecheck_select(mrna, agn_typecheck_exon);
  if(gt_array_size(s->exons) > 1)
    gt_array_sort(s->exons, (GtCompare)agn_genome_node_compare);
  s->exon_tree = gt_interval_tree_new(NULL);
  for(i = 0; i < gt_array_size(s->exons); i++)
  {
    GtFeatureNode *exon = *(GtFeatureNode **)gt_array_get(s->exons, i);
    GtRange r = gt_genome_node_get_range((GtGenomeNode *)exon);
    GtIntervalTreeNode *itn = gt_interval_tree_node_new(exon, r.start, r.end);
    gt_interval_tree_insert(s->exon_tree, itn);
  }
  return s;
}

static void mrna_structure_delete(MrnaStructure *s)
{
  gt_interval_tree_delete(s->exon_tree);
  gt_free(s);
}

static void check_alt_splice_sites(MrnaStructure *s1, MrnaStructure *s2)
{
  GtUword i;
  if(gt_array_size(s1->exons) <= 1 || gt_array_size(s2->exons) <= 1)
    return;

  // Left to right
  for(i = 1; i < gt_array_size(s1->exons); i++)
  {
    GtFeatureNode *exon_1a = *(GtFeatureNode **)gt_array_get(s1->exons, i-1);
    GtRange r1a = gt_genome_node_get_range((GtGenomeNode *)exon_1a);
    GtFeatureNode *exon_1b = *(GtFeatureNode **)gt_array_get(s1->exons, i);
    GtRange r1b = gt_genome_node_get_range((GtGenomeNode *)exon_1b);

    GtIntervalTreeNode *exon_2a_itn = gt_interval_tree_find_first_overlapping(
                                          s2->exon_tree, r1a.start, r1a.end);
    GtIntervalTreeNode *exon_2b_itn = gt_interval_tree_find_first_overlapping(
                                          s2->exon_tree, r1a.start, r1a.end);
    if(exon_2a_itn == NULL || exon_2b_itn == NULL)
      continue;
    GtFeatureNode *exon_2a = gt_interval_tree_node_get_data(exon_2a_itn);
    GtRange r2a = gt_genome_node_get_range((GtGenomeNode *)exon_2a);
    GtFeatureNode *exon_2b = gt_interval_tree_node_get_data(exon_2b_itn);
    GtRange r2b = gt_genome_node_get_range((GtGenomeNode *)exon_2b);

    if(r1a.start == r2a.start && r1a.end != r2a.end && r1b.start == r2b.start)
    {
      GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)exon_1a);
      GtStrand strand = gt_feature_node_get_strand(exon_1a);
      const char *sitetype = "donor";
      if(strand == GT_STRAND_REVERSE)
        sitetype = "acceptor";
      printf("alternate %s site: %s [%lu, %lu] vs [%lu, %lu]\n", sitetype,
             gt_str_get(seqid), r1a.start, r1a.end, r2a.start, r2a.end);
    }
  }

  // Right to left
  for(i = gt_array_size(s1->exons) - 1; i > 0; i--)
  {
    GtFeatureNode *exon_1a = *(GtFeatureNode **)gt_array_get(s1->exons, i);
    GtRange r1a = gt_genome_node_get_range((GtGenomeNode *)exon_1a);
    GtFeatureNode *exon_1b = *(GtFeatureNode **)gt_array_get(s1->exons, i-1);
    GtRange r1b = gt_genome_node_get_range((GtGenomeNode *)exon_1b);

    GtIntervalTreeNode *exon_2a_itn = gt_interval_tree_find_first_overlapping(
                                          s2->exon_tree, r1a.start, r1a.end);
    GtIntervalTreeNode *exon_2b_itn = gt_interval_tree_find_first_overlapping(
                                          s2->exon_tree, r1a.start, r1a.end);
    if(exon_2a_itn == NULL || exon_2b_itn == NULL)
      continue;
    GtFeatureNode *exon_2a = gt_interval_tree_node_get_data(exon_2a_itn);
    GtRange r2a = gt_genome_node_get_range((GtGenomeNode *)exon_2a);
    GtFeatureNode *exon_2b = gt_interval_tree_node_get_data(exon_2b_itn);
    GtRange r2b = gt_genome_node_get_range((GtGenomeNode *)exon_2b);

    if(r1a.end == r2a.end && r1a.start != r2a.start && r1b.end == r2b.end)
    {
      GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)exon_1a);
      GtStrand strand = gt_feature_node_get_strand(exon_1a);
      const char *sitetype = "acceptor";
      if(strand == GT_STRAND_REVERSE)
        sitetype = "donor";
      printf("alternate %s site: %s [%lu, %lu] vs [%lu, %lu]\n", sitetype,
             gt_str_get(seqid), r1a.start, r1a.end, r2a.start, r2a.end);
    }
  }
}

static void check_skipped_exons(MrnaStructure *s1, MrnaStructure *s2)
{
  GtUword i;
  if(gt_array_size(s1->exons) <= 1 || gt_array_size(s2->exons) <= 1)
    return;

  for(i = 1; i < gt_array_size(s1->exons); i++)
  {
    GtFeatureNode *exon_1a = *(GtFeatureNode **)gt_array_get(s1->exons, i-1);
    GtRange r1a = gt_genome_node_get_range((GtGenomeNode *)exon_1a);
    GtFeatureNode *exon_1b = *(GtFeatureNode **)gt_array_get(s1->exons, i);
    GtRange r1b = gt_genome_node_get_range((GtGenomeNode *)exon_1b);

    GtArray *exons2 = gt_array_new(sizeof(GtIntervalTreeNode *));
    gt_interval_tree_find_all_overlapping(s2->exon_tree, r1a.start, r1b.end,
                                          exons2);
    if(gt_array_size(exons2) > 2)
    {
      gt_array_sort(exons2, (GtCompare)agn_genome_node_compare);
      GtFeatureNode *exon_2a = *(GtFeatureNode **)gt_array_get_first(exons2);
      GtFeatureNode *exon_2b = *(GtFeatureNode **)gt_array_get_last(exons2);
      GtRange r2a = gt_genome_node_get_range((GtGenomeNode *)exon_2a);
      GtRange r2b = gt_genome_node_get_range((GtGenomeNode *)exon_2b);
      if(gt_range_compare(&r1a, &r2a) == 0 && gt_range_compare(&r1b, &r2b) == 0)
      {
        GtUword numskipped = gt_array_size(exons2) - 2;
        GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)exon_1a);
        printf("%lu skipped exon(s): %s [%lu, %lu]\n", numskipped,
               gt_str_get(seqid), r1a.start, r1b.end);
      }
    }
    gt_array_delete(exons2);
  }

  for(i = 1; i < gt_array_size(s2->exons); i++)
  {
    GtFeatureNode *exon_1a = *(GtFeatureNode **)gt_array_get(s2->exons, i-1);
    GtRange r1a = gt_genome_node_get_range((GtGenomeNode *)exon_1a);
    GtFeatureNode *exon_1b = *(GtFeatureNode **)gt_array_get(s2->exons, i);
    GtRange r1b = gt_genome_node_get_range((GtGenomeNode *)exon_1b);

    GtArray *exons2 = gt_array_new(sizeof(GtIntervalTreeNode *));
    gt_interval_tree_find_all_overlapping(s1->exon_tree, r1a.start, r1b.end,
                                          exons2);
    if(gt_array_size(exons2) > 2)
    {
      gt_array_sort(exons2, (GtCompare)agn_genome_node_compare);
      GtFeatureNode *exon_2a = *(GtFeatureNode **)gt_array_get_first(exons2);
      GtFeatureNode *exon_2b = *(GtFeatureNode **)gt_array_get_last(exons2);
      GtRange r2a = gt_genome_node_get_range((GtGenomeNode *)exon_2a);
      GtRange r2b = gt_genome_node_get_range((GtGenomeNode *)exon_2b);
      if(gt_range_compare(&r1a, &r2a) == 0 && gt_range_compare(&r1b, &r2b) == 0)
      {
        GtUword numskipped = gt_array_size(exons2) - 2;
        GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)exon_1a);
        printf("%lu skipped exon(s): %s [%lu, %lu]\n", numskipped,
               gt_str_get(seqid), r1a.start, r1b.end);
      }
    }
    gt_array_delete(exons2);
  }
}

int main(int argc, char **argv)
{
  GtError *error;
  GtGenomeNode *gn;
  GtLogger *logger;
  GtNodeStream *current_stream, *previous_stream;
  GtQueue *streams;
  int had_err;

  gt_lib_init();
  streams = gt_queue_new();

  current_stream = gt_gtf_in_stream_new(argv[1]);
  gt_queue_add(streams, current_stream);
  previous_stream = current_stream;

  logger = gt_logger_new(true, "", stderr);
  current_stream = agn_locus_stream_new(previous_stream, logger);
  gt_queue_add(streams, current_stream);
  previous_stream = current_stream;

  error = gt_error_new();
  while (!(had_err = gt_node_stream_next(previous_stream, &gn, error)) && gn)
  {
    GtFeatureNode *fn = gt_feature_node_try_cast(gn);
    if(fn && gt_feature_node_has_type(fn, "locus"))
    {
      GtUword i,j, numrnas;
      GtArray *mrnas = agn_typecheck_select(fn, agn_typecheck_mrna);
      numrnas = gt_array_size(mrnas);
      if(numrnas <= 1)
      {
        gt_array_delete(mrnas);
        gt_genome_node_delete(gn);
        continue;
      }
      for(i = 0; i < numrnas; i++)
      {
        GtFeatureNode *mrnai = *(GtFeatureNode **)gt_array_get(mrnas, i);
        for(j = i+1; j < numrnas; j++)
        {
          GtFeatureNode *mrnaj = *(GtFeatureNode **)gt_array_get(mrnas, j);

          MrnaStructure *s_i = mrna_structure_new(mrnai);
          MrnaStructure *s_j = mrna_structure_new(mrnaj);
          check_alt_splice_sites(s_i, s_j);
          check_skipped_exons(s_i, s_j);
          mrna_structure_delete(s_i);
          mrna_structure_delete(s_j);
        }
      }
    }
    gt_genome_node_delete(gn);
  }
  if(had_err)
  {
    fprintf(stderr, "error: %s\n", gt_error_get(error));
  }

  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  gt_error_delete(error);
  gt_logger_delete(logger);
  gt_lib_clean();
  return 0;
}
