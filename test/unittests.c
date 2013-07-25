#include <string.h>
#include "AgnGeneLocus.h"
#include "AgnUtils.h"
#include "AgnTranscriptClique.h"

GtFeatureNode *agn_unit_test_eden()
{
  GtGenomeNode *gene, *mrna1, *mrna2, *mrna3, *feature;
  GtGenomeNode *f1, *f2, *f3, *f4;
  GtStr *seqid = gt_str_new_cstr("ctg123");

  gene = gt_feature_node_new(seqid, "gene", 1000, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_attribute((GtFeatureNode *)gene, "ID", "EDEN");

  mrna1 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_attribute((GtFeatureNode *)mrna1, "ID", "EDEN.1");
  gt_feature_node_add_child((GtFeatureNode *)gene, (GtFeatureNode *)mrna1);
  mrna2 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_attribute((GtFeatureNode *)mrna2, "ID", "EDEN.2");
  gt_feature_node_add_child((GtFeatureNode *)gene, (GtFeatureNode *)mrna2);
  mrna3 = gt_feature_node_new(seqid, "mRNA", 1300, 9000, GT_STRAND_FORWARD);
  gt_feature_node_add_attribute((GtFeatureNode *)mrna3, "ID", "EDEN.3");
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
  f4 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
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
  puts("    AgnTranscriptClique class");

  GtFeatureNode *eden = agn_unit_test_eden();
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)eden);
  AgnGeneLocus *locus = agn_gene_locus_new(gt_str_get(seqid));
  agn_gene_locus_add_gene(locus, eden);

  GtArray *trans = agn_gene_locus_get_transcripts(locus);
  GtArray *cliques = agn_enumerate_feature_cliques(trans);

  bool parsearraypass = gt_array_size(cliques) == 3;
  if(parsearraypass)
    printf("        | %-30s| PASS\n", "parse from array");
  else
    printf("        | %-30s| FAIL\n", "parse from array");
  
  bool numtranspass = true;
  bool cdslenpass =  true;
  bool iterpass = true;
  unsigned long i;
  for(i = 0; i < gt_array_size(cliques); i++)
  {
    AgnTranscriptClique *tc = *(AgnTranscriptClique **)gt_array_get(cliques, i);
    if(agn_transcript_clique_size(tc) != 1)
    {
      numtranspass = false;
    }
    GtFeatureNode *mrna = agn_transcript_clique_next(tc);
    const char *mrnaid = gt_feature_node_get_attribute(mrna, "ID");
    unsigned long cdslength = agn_transcript_clique_cds_length(tc);
    if(strcmp(mrnaid, "EDEN.1") == 0)
    {
      if(cdslength != 2305) cdslenpass = false;
    }
    else if(strcmp(mrnaid, "EDEN.2") == 0)
    {
      if(cdslength != 1402) cdslenpass = false;
    }
    else if(strcmp(mrnaid, "EDEN.3") == 0)
    {
      if(cdslength != 1704) cdslenpass = false;
    }
    else
    {
      iterpass = false;
    }
  }
  if(numtranspass)
    printf("        | %-30s| PASS\n", "transcripts per clique");
  else
    printf("        | %-30s| FAIL\n", "transcripts per clique");
  if(cdslenpass)
    printf("        | %-30s| PASS\n", "clique CDS length");
  else
    printf("        | %-30s| FAIL\n", "clique CDS length");
  if(iterpass)
    printf("        | %-30s| PASS\n", "iterate through transcripts");
  else
    printf("        | %-30s| FAIL\n", "iterate through transcripts");

  return parsearraypass && numtranspass && cdslenpass && iterpass;
}

int main(int argc, char **argv)
{
  puts("AEGeAn Unit Tests");
  gt_lib_init();
  bool tc_pass = agn_transcript_clique_unit_test();
  if(tc_pass){ /* this line for compiler warning */ }
  gt_lib_clean();
  return 0;
}
