#include "AgnGeneLocus.h"
#include "AgnUtils.h"
#include "AgnTranscriptClique.h"

GtFeatureNode *agn_unit_test_eden()
{
  GtGenomeNode *gene, *mrna1, *mrna2, *mrna3, *feature;
  GtGenomeNode *f1, *f2, *f3, *f4;
  GtStr *seqid = gt_str_new_cstr("ctg123");

  gene = gt_feature_node_new(seqid, "gene", 1000, 9000, GT_STRAND_FORWARD);

  mrna1 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene, (GtFeatureNode *)mrna1);
  mrna2 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene, (GtFeatureNode *)mrna2);
  mrna3 = gt_feature_node_new(seqid, "mRNA", 1300, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene, (GtFeatureNode *)mrna3);

  feature = gt_feature_node_new(seqid, "exon", 1050, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)feature);
  feature = gt_feature_node_new(seqid, "exon", 1300, 1500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)feature);
  feature = gt_feature_node_new(seqid, "exon", 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)feature);
  feature = gt_feature_node_new(seqid, "exon", 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)feature);
  feature = gt_feature_node_new(seqid, "exon", 7000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)feature);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)feature);

  f1 = gt_feature_node_new(seqid, "CDS", 1201, 1500, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f1);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f1, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)f1);
  f2 = gt_feature_node_new(seqid, "CDS", 3000, 3902, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f2);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f2, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)f2);
  f3 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f3);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f3, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)f3);
  f4 = gt_feature_node_new(seqid, "CDS", 7000, 76000, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f4);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f4, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)f4);

  f1 = gt_feature_node_new(seqid, "CDS", 1201, 1500, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f1);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f1, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)f1);
  f2 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f2);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f2, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)f2);
  f3 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f3);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f3, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)f3);

  f1 = gt_feature_node_new(seqid, "CDS", 3301, 3902, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f1);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f1, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)f1);
  f2 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f2);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f2, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)f2);
  f3 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
  gt_feature_node_make_multi_representative((GtFeatureNode *)f3);
  gt_feature_node_set_multi_representative((GtFeatureNode *)f3, (GtFeatureNode *)f1);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)f3);

  return (GtFeatureNode *)gene;
}

bool agn_transcript_clique_unit_test()
{
  puts("[AgnTranscriptClique::UnitTest]");

  GtFeatureNode *eden = agn_unit_test_eden();
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)eden);
  AgnGeneLocus *locus = agn_gene_locus_new(gt_str_get(seqid));
  agn_gene_locus_add_gene(locus, eden);

  GtArray *trans = agn_gene_locus_get_transcripts(locus);
  GtArray *cliques = agn_enumerate_feature_cliques(trans);

  bool parsearraypass = gt_array_size(cliques) == 3;
  if(parsearraypass)
    printf("    | %-25s| PASS\n", "parse from array");
  else
    printf("    | %-25s| FAIL\n", "parse from array");
  
  bool numtranspass = true;
  unsigned long i;
  for(i = 0; i < gt_array_size(cliques); i++)
  {
    AgnTranscriptClique *tc = *(AgnTranscriptClique **)gt_array_get(cliques, i);
    if(agn_transcript_clique_size(tc) != 1)
    {
      numtranspass = false;
      break;
    }
  }
  if(numtranspass)
    printf("    | %-25s| PASS\n", "transcripts per clique");
  else
    printf("    | %-25s| FAIL\n", "transcripts per clique");

  return parsearraypass && numtranspass;
}

int main(int argc, char **argv)
{
  gt_lib_init();
  bool tc_pass = agn_transcript_clique_unit_test();
  if(tc_pass){ /* this line for compiler warning */ }
  gt_lib_clean();
  return 0;
}
