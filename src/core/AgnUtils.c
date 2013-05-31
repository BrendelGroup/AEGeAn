#include <string.h>
#include <time.h>
#include "extended/feature_node.h"
#include "AgnCanonNodeVisitor.h"
#include "AgnGeneValidator.h"
#include "AgnGtExtensions.h"
#include "AgnPairwiseCompareLocus.h"
#include "AgnUtils.h"

void agn_bron_kerbosch( GtArray *R, GtArray *P, GtArray *X, GtArray *cliques,
                        bool skipsimplecliques )
{
  gt_assert(R != NULL && P != NULL && X != NULL && cliques != NULL);

  if(gt_array_size(P) == 0 && gt_array_size(X) == 0)
  {
    if(skipsimplecliques == false || gt_array_size(R) != 1)
    {
      unsigned long i;
      AgnTranscriptClique *clique = agn_transcript_clique_new();
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
    GtArray *newR = agn_gt_array_copy(R, sizeof(GtGenomeNode *));
    gt_array_add(newR, v);
    // newP = P \intersect N(v)
    GtArray *newP = agn_feature_neighbors(v, P);
    // newX = X \intersect N(v)
    GtArray *newX = agn_feature_neighbors(v, X);

    // Recursive call
    // agn_bron_kerbosch(R \union {v}, P \intersect N(v), X \intersect N(X))
    agn_bron_kerbosch(newR, newP, newX, cliques, skipsimplecliques);

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

double agn_calc_edit_distance(GtFeatureNode *t1, GtFeatureNode *t2)
{
// fprintf(stderr, "DELETEME newtest a\n");
  AgnTranscriptClique *clique1 = agn_transcript_clique_new();
  agn_transcript_clique_add(clique1, t1);
  AgnTranscriptClique *clique2 = agn_transcript_clique_new();
  agn_transcript_clique_add(clique2, t2);
// fprintf(stderr, "DELETEME newtest b\n");

  GtGenomeNode *gn1 = (GtGenomeNode *)t1;
  GtRange r1 = gt_genome_node_get_range(gn1);
  GtStr *seqid = gt_genome_node_get_seqid(gn1);
  GtGenomeNode *gn2 = (GtGenomeNode *)t2;
  GtRange r2 = gt_genome_node_get_range(gn2);
  GtRange local_range = r1;
  if(r2.start < r1.start)
    local_range.start = r2.start;
  if(r2.end > r1.end)
    local_range.end = r2.end;
// fprintf(stderr, "DELETEME newtest c\n");

  AgnCliquePair *pair = agn_clique_pair_new(gt_str_get(seqid), clique1, clique2,
                                            &local_range);
  agn_clique_pair_build_model_vectors(pair);
  agn_clique_pair_comparative_analysis(pair);
// fprintf(stderr, "DELETEME newtest d\n");

  double ed = agn_clique_pair_get_edit_distance(pair);
// fprintf(stderr, "DELETEME newtest e\n");

  agn_transcript_clique_delete(clique1);
// fprintf(stderr, "DELETEME newtest f\n");
  agn_transcript_clique_delete(clique2);
// fprintf(stderr, "DELETEME newtest g\n");
  agn_clique_pair_delete(pair);
// fprintf(stderr, "DELETEME newtest h\n");

  return ed;
}

double agn_calc_splice_complexity(GtArray *transcripts)
{
  unsigned long n = gt_array_size(transcripts);
  unsigned long i,j;
  double sc = 0.0;

  for(i = 0; i < n; i++)
  {
    GtFeatureNode *t_i = *(GtFeatureNode **)gt_array_get(transcripts, i);
    for(j = 0; j < i; j++)
    {
      GtFeatureNode *t_j = *(GtFeatureNode **)gt_array_get(transcripts, j);
// fprintf(stderr, "DELETEME t_i,t_j=%p,%p\n", t_i, t_j);
      if(agn_gt_feature_node_overlap(t_i, t_j))
      {
// fprintf(stderr, "DELETEME i,j=%lu,%lu\n", i, j);
        sc += agn_calc_edit_distance(t_i, t_j);
      }
    }
  }

  return sc;
}

GtArray* agn_enumerate_feature_cliques(GtArray *feature_set)
{
  GtArray *cliques = gt_array_new( sizeof(GtArray *) );
  
  if(gt_array_size(feature_set) == 1)
  {
    GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(feature_set, 0);
    AgnTranscriptClique *clique = agn_transcript_clique_new();
    agn_transcript_clique_add(clique, fn);
    gt_array_add(cliques, clique);
  }
  else
  {
    // First add each transcript as a clique, even if it is not a maximal clique
    unsigned long i;
    for(i = 0; i < gt_array_size(feature_set); i++)
    {
      GtFeatureNode *fn = *(GtFeatureNode **)gt_array_get(feature_set, i);
      AgnTranscriptClique *clique = agn_transcript_clique_new();
      agn_transcript_clique_add(clique, fn);
      gt_array_add(cliques, clique);
    }
    
    // Then use the Bron-Kerbosch algorithm to find all maximal cliques
    // containing >1 transcript
    GtArray *R = gt_array_new( sizeof(GtGenomeNode *) );
    GtArray *P = agn_gt_array_copy(feature_set, sizeof(GtGenomeNode *));
    GtArray *X = gt_array_new( sizeof(GtGenomeNode *) );
  
    // Initial call: agn_bron_kerbosch(\emptyset, vertex_set, \emptyset )
    agn_bron_kerbosch(R, P, X, cliques, true);
  
    gt_array_delete(R);
    gt_array_delete(P);
    gt_array_delete(X);
  }

  return cliques;
}

GtArray* agn_feature_neighbors(GtGenomeNode *feature, GtArray *feature_set)
{
  GtArray *neighbors = gt_array_new( sizeof(GtGenomeNode *) );
  unsigned long i;
  for(i = 0; i < gt_array_size(feature_set); i++)
  {
    GtGenomeNode *other = *(GtGenomeNode **)gt_array_get(feature_set, i);
    if(other != feature)
    {
      GtRange feature_range = gt_genome_node_get_range(feature);
      GtRange other_range = gt_genome_node_get_range(other);
      if(gt_range_overlap(&feature_range, &other_range) == false)
        gt_array_add(neighbors, other);
    }
  }
  return neighbors;
}

FILE *agn_fopen(const char *filename, const char *mode)
{
  FILE *fp = fopen(filename, mode);
  if(fp == NULL)
  {
    fprintf(stderr, "error: could not open '%s'\n", filename);
    exit(1);
  }
  return fp;
}

GtFeatureIndex *agn_import_canonical(int numfiles, const char **filenames,
                                     AgnLogger *logger)
{
  GtNodeStream *gff3 = gt_gff3_in_stream_new_unsorted(numfiles, filenames);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3);
  
  GtFeatureIndex *features = gt_feature_index_memory_new();
  AgnGeneValidator *validator = agn_gene_validator_new();
  GtNodeVisitor *nv = agn_canon_node_visitor_new(features, validator, logger);

  GtGenomeNode *gn;
  bool loaderror;
  GtError *error = gt_error_new();
  while(!(loaderror = gt_node_stream_next(gff3, &gn, error)) && gn)
  {
    gt_genome_node_accept(gn, nv, error);
    if(agn_logger_has_error(logger))
      break;
  }
  gt_node_stream_delete(gff3);
  gt_node_visitor_delete(nv);
  agn_gene_validator_delete(validator);

  if(agn_logger_has_error(logger))
  {
    gt_feature_index_delete(features);
    features = NULL;
  }
  
  gt_error_delete(error);
  return features;
}

bool agn_infer_cds_range_from_exon_and_codons(GtRange *exon_range,
                                              GtRange *leftcodon_range,
                                              GtRange *rightcodon_range,
                                              GtRange *cds_range)
{
  cds_range->start = 0;
  cds_range->end   = 0;
  
  // UTR
  if(exon_range->end < leftcodon_range->start ||
     exon_range->start > rightcodon_range->end)
    return false;
  
  bool overlap_left  = gt_range_overlap(exon_range, leftcodon_range);
  bool overlap_right = gt_range_overlap(exon_range, rightcodon_range);
  if(overlap_left && overlap_right)
  {
    cds_range->start = leftcodon_range->start;
    cds_range->end   = rightcodon_range->end;
  }
  else if(overlap_left)
  {
    cds_range->start = leftcodon_range->start;
    cds_range->end   = exon_range->end;
  }
  else if(overlap_right)
  {
    cds_range->start = exon_range->start;
    cds_range->end   = rightcodon_range->end;
  }
  else
  {
    cds_range->start = exon_range->start;
    cds_range->end   = exon_range->end;
  }
  
  return true;
}

GtStrArray* agn_seq_intersection(GtFeatureIndex *refrfeats,
                                 GtFeatureIndex *predfeats, AgnLogger *logger)
{
  // Fetch seqids from reference and prediction annotations
  GtError *e = gt_error_new();
  GtStrArray *refrseqids = gt_feature_index_get_seqids(refrfeats, e);
  if(gt_error_is_set(e))
  {
    agn_logger_log_error(logger, "error fetching seqids for reference: %s",
                         gt_error_get(e));
    gt_error_unset(e);
  }
  GtStrArray *predseqids = gt_feature_index_get_seqids(predfeats, e);
  if(gt_error_is_set(e))
  {
    agn_logger_log_error(logger, "error fetching seqids for prediction: %s",
                         gt_error_get(e));
    gt_error_unset(e);
  }
  gt_error_delete(e);
  if(agn_logger_has_error(logger))
  {
    gt_str_array_delete(refrseqids);
    gt_str_array_delete(predseqids);
    return NULL;
  }
  GtStrArray *seqids = agn_gt_str_array_intersection(refrseqids, predseqids);
  
  // Print reference sequences with no prediction annotations
  unsigned long i, j;
  for(i = 0; i < gt_str_array_size(refrseqids); i++)
  {
    const char *refrseq = gt_str_array_get(refrseqids, i);
    int matches = 0;
    for(j = 0; j < gt_str_array_size(seqids); j++)
    {
      const char *seq = gt_str_array_get(seqids, j);
      if(strcmp(refrseq, seq) == 0)
        matches++;
    }
    if(matches == 0)
    {
      agn_logger_log_warning(logger, "no prediction annotations found for "
                             "sequence '%s'", refrseq);
    }
  }

  // Print prediction sequences with no reference annotations
  for(i = 0; i < gt_str_array_size(predseqids); i++)
  {
    const char *predseq = gt_str_array_get(predseqids, i);
    int matches = 0;
    for(j = 0; j < gt_str_array_size(seqids); j++)
    {
      const char *seq = gt_str_array_get(seqids, j);
      if(strcmp(predseq, seq) == 0)
        matches++;
    }
    if(matches == 0)
    {
      agn_logger_log_warning(logger, "no reference annotations found for "
                             "sequence '%s'", predseq);
    }
  }
  
  if(gt_str_array_size(seqids) == 0)
  {
    agn_logger_log_error(logger, "no sequences in common between reference and "
                         "prediction");
  }

  gt_str_array_delete(refrseqids);
  gt_str_array_delete(predseqids);
  return seqids;
}

GtStrArray* agn_seq_union(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats,
                          AgnLogger *logger)
{
  // Fetch seqids from reference and prediction annotations
  GtError *e = gt_error_new();
  GtStrArray *refrseqids = gt_feature_index_get_seqids(refrfeats, e);
  if(gt_error_is_set(e))
  {
    agn_logger_log_error(logger, "error fetching seqids for reference: %s",
                         gt_error_get(e));
    gt_error_unset(e);
  }
  GtStrArray *predseqids = gt_feature_index_get_seqids(predfeats, e);
  if(gt_error_is_set(e))
  {
    agn_logger_log_error(logger, "error fetching seqids for prediction: %s",
                         gt_error_get(e));
    gt_error_unset(e);
  }
  gt_error_delete(e);
  if(agn_logger_has_error(logger))
  {
    gt_str_array_delete(refrseqids);
    gt_str_array_delete(predseqids);
    return NULL;
  }
  GtStrArray *seqids = agn_gt_str_array_union(refrseqids, predseqids);
  
  gt_str_array_delete(refrseqids);
  gt_str_array_delete(predseqids);
  return seqids;
}

int agn_sprintf_comma(unsigned long n, char *buffer)
{
  if(n < 1000)
  {
    int spaces = sprintf(buffer, "%lu", n);
    buffer += spaces;
    return spaces;
  }
  int spaces = agn_sprintf_comma(n / 1000, buffer);
  buffer += spaces;
  sprintf(buffer, ",%03lu", n % 1000);
  return spaces + 4;
}

int agn_string_compare(const void *p1, const void *p2)
{
  const char *s1 = *(char **)p1;
  const char *s2 = *(char **)p2;
  return strcmp(s1, s2);
}

GtRange agn_transcript_cds_range(GtFeatureNode *transcript)
{
  gt_assert(transcript);
  GtRange trange;
  trange.start = 0;
  trange.end = 0;
  
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(transcript);
  GtFeatureNode *current;
  for
  (
    current = gt_feature_node_iterator_next(iter);
    current != NULL;
    current = gt_feature_node_iterator_next(iter)
  )
  {
    if(agn_gt_feature_node_is_cds_feature(current))
    {
      GtRange crange = gt_genome_node_get_range((GtGenomeNode *)current);
      if(trange.start == 0 || crange.start < trange.start)
        trange.start = crange.start;
      if(trange.end == 0 || crange.end > trange.end)
        trange.end = crange.end;
    }
  }
  
  if(gt_feature_node_get_strand(transcript) == GT_STRAND_REVERSE)
  {
    unsigned long temp = trange.start;
    trange.start = trange.end;
    trange.end = temp;
  }
  return trange;
}

void agn_transcript_structure_gbk(GtFeatureNode *transcript, FILE *outstream)
{
  gt_assert(transcript && outstream);
  
  GtArray *exons = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(transcript);
  GtFeatureNode *child;
  for
  (
    child = gt_feature_node_iterator_next(iter);
    child != NULL;
    child = gt_feature_node_iterator_next(iter)
  )
  {
    if(agn_gt_feature_node_is_exon_feature(child))
      gt_array_add(exons, child);
  }
  gt_feature_node_iterator_delete(iter);
  
  gt_assert(gt_array_size(exons) > 0);
  gt_array_sort(exons, (GtCompare)agn_gt_genome_node_compare);
  
  if(gt_feature_node_get_strand(transcript) == GT_STRAND_REVERSE)
    fputs("complement(", outstream);
  
  if(gt_array_size(exons) == 1)
  {
    GtGenomeNode *exon = *(GtGenomeNode **)gt_array_get(exons, 0);
    GtRange exonrange = gt_genome_node_get_range(exon);
    fprintf(outstream, "<%lu..>%lu", exonrange.start, exonrange.end);
  }
  else
  {
    fputs("join(", outstream);
    unsigned long i;
    for(i = 0; i < gt_array_size(exons); i++)
    {
      GtGenomeNode *exon = *(GtGenomeNode **)gt_array_get(exons, i);
      GtRange exonrange = gt_genome_node_get_range(exon);
      
      if(i == 0)
        fprintf(outstream, "<%lu..%lu", exonrange.start, exonrange.end);
      else if(i+1 == gt_array_size(exons))
        fprintf(outstream, ",%lu..>%lu", exonrange.start, exonrange.end);
      else
        fprintf(outstream, ",%lu..%lu", exonrange.start, exonrange.end);
    }
    fputs(")", outstream);
  }
  
  if(gt_feature_node_get_strand(transcript) == GT_STRAND_REVERSE)
    fputs(")", outstream);
}
