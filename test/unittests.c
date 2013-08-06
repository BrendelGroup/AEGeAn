#include <string.h>
#include "AgnGeneLocus.h"
#include "AgnUtils.h"
#include "AgnTranscriptClique.h"

GtFeatureNode *agn_unit_test_eden()
{
  GtGenomeNode *gene, *mrna1, *mrna2, *mrna3, *feature;
  GtFeatureNode *genefn, *mrna1fn, *mrna2fn, *mrna3fn, *featurefn;
  GtGenomeNode *f1, *f2, *f3, *f4;
  GtFeatureNode *fn1, *fn2, *fn3, *fn4;
  GtStr *seqid = gt_str_new_cstr("ctg123");

  gene = gt_feature_node_new(seqid, "gene", 1000, 9000, GT_STRAND_FORWARD);
  genefn = (GtFeatureNode *)gene;
  gt_feature_node_add_attribute(genefn, "ID", "EDEN");

  mrna1 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  mrna1fn = (GtFeatureNode *)mrna1;
  gt_feature_node_add_attribute(mrna1fn, "ID", "EDEN.1");
  gt_feature_node_add_child(genefn, mrna1fn);
  mrna2 = gt_feature_node_new(seqid, "mRNA", 1050, 9000, GT_STRAND_FORWARD);
  mrna2fn = (GtFeatureNode *)mrna2;
  gt_feature_node_add_attribute(mrna2fn, "ID", "EDEN.2");
  gt_feature_node_add_child(genefn, mrna2fn);
  mrna3 = gt_feature_node_new(seqid, "mRNA", 1300, 9000, GT_STRAND_FORWARD);
  mrna3fn = (GtFeatureNode *)mrna3;
  gt_feature_node_add_attribute(mrna3fn, "ID", "EDEN.3");
  gt_feature_node_add_child(genefn, mrna3fn);

  feature = gt_feature_node_new(seqid, "exon", 1050, 1500, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna2fn, featurefn);
  feature = gt_feature_node_new(seqid, "exon", 1300, 1500, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna3fn, featurefn);
  feature = gt_feature_node_new(seqid, "exon", 3000, 3902, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);
  feature = gt_feature_node_new(seqid, "exon", 5000, 5500, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna2fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);
  feature = gt_feature_node_new(seqid, "exon", 7000, 9000, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna2fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);

  f1 = gt_feature_node_new(seqid, "CDS", 1201, 1500, GT_STRAND_FORWARD);
  fn1 = (GtFeatureNode *)f1;
  gt_feature_node_make_multi_representative(fn1);
  gt_feature_node_set_multi_representative(fn1, fn1);
  gt_feature_node_add_child(mrna1fn, fn1);
  f2 = gt_feature_node_new(seqid, "CDS", 3000, 3902, GT_STRAND_FORWARD);
  fn2 = (GtFeatureNode *)f2;
  gt_feature_node_make_multi_representative(fn2);
  gt_feature_node_set_multi_representative(fn2, fn1);
  gt_feature_node_add_child(mrna1fn, fn2);
  f3 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  fn3 = (GtFeatureNode *)f3;
  gt_feature_node_make_multi_representative(fn3);
  gt_feature_node_set_multi_representative(fn3, fn1);
  gt_feature_node_add_child(mrna1fn, fn3);
  f4 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
  fn4 = (GtFeatureNode *)f4;
  gt_feature_node_make_multi_representative(fn4);
  gt_feature_node_set_multi_representative(fn4, fn1);
  gt_feature_node_add_child(mrna1fn, fn4);

  f1 = gt_feature_node_new(seqid, "CDS", 1201, 1500, GT_STRAND_FORWARD);
  fn1 = (GtFeatureNode *)f1;
  gt_feature_node_make_multi_representative(fn1);
  gt_feature_node_set_multi_representative(fn1, fn1);
  gt_feature_node_add_child(mrna2fn, fn1);
  f2 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  fn2 = (GtFeatureNode *)f2;
  gt_feature_node_make_multi_representative(fn2);
  gt_feature_node_set_multi_representative(fn2, fn1);
  gt_feature_node_add_child(mrna2fn, fn2);
  f3 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
  fn3 = (GtFeatureNode *)f3;
  gt_feature_node_make_multi_representative(fn3);
  gt_feature_node_set_multi_representative(fn3, fn1);
  gt_feature_node_add_child(mrna2fn, fn3);

  f1 = gt_feature_node_new(seqid, "CDS", 3301, 3902, GT_STRAND_FORWARD);
  fn1 = (GtFeatureNode *)f1;
  gt_feature_node_make_multi_representative(fn1);
  gt_feature_node_set_multi_representative(fn1, fn1);
  gt_feature_node_add_child(mrna3fn, fn1);
  f2 = gt_feature_node_new(seqid, "CDS", 5000, 5500, GT_STRAND_FORWARD);
  fn2 = (GtFeatureNode *)f2;
  gt_feature_node_make_multi_representative(fn2);
  gt_feature_node_set_multi_representative(fn2, fn1);
  gt_feature_node_add_child(mrna3fn, fn2);
  f3 = gt_feature_node_new(seqid, "CDS", 7000, 7600, GT_STRAND_FORWARD);
  fn3 = (GtFeatureNode *)f3;
  gt_feature_node_make_multi_representative(fn3);
  gt_feature_node_set_multi_representative(fn3, fn1);
  gt_feature_node_add_child(mrna3fn, fn3);

  return genefn;
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
    const char *mrnaid = agn_transcript_clique_id(tc);
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
