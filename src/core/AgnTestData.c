#include "AgnTestData.h"

GtArray *agn_test_data_genes_codons()
{
  GtArray *genes = gt_array_new( sizeof(GtGenomeNode *) );
  GtStr *seqid = gt_str_new_cstr("chr8");

  GtGenomeNode *gene1  = gt_feature_node_new(seqid, "gene", 22057, 23119,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *mrna1  = gt_feature_node_new(seqid, "mRNA", 22057, 23119,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *exon1  = gt_feature_node_new(seqid, "exon", 22057, 22382,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *exon2  = gt_feature_node_new(seqid, "exon", 22497, 22550,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *exon3  = gt_feature_node_new(seqid, "exon", 22651, 23119,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *start1 = gt_feature_node_new(seqid, "start_codon", 22167, 22169,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *stop1  = gt_feature_node_new(seqid, "stop_codon", 23020, 23022,
                                             GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene1, (GtFeatureNode *)mrna1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)exon1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)exon2);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)exon3);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)start1);
  gt_feature_node_add_child((GtFeatureNode *)mrna1, (GtFeatureNode *)stop1);
  gt_array_add(genes, gene1);

  GtGenomeNode *gene2  = gt_feature_node_new(seqid, "gene", 48012, 48984,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *mrna2  = gt_feature_node_new(seqid, "mRNA", 48012, 48984,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon4  = gt_feature_node_new(seqid, "exon", 48012, 48537,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon5  = gt_feature_node_new(seqid, "exon", 48637, 48766,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon6  = gt_feature_node_new(seqid, "exon", 48870, 48984,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *start2 = gt_feature_node_new(seqid, "start_codon", 48982, 48984,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *stop2  = gt_feature_node_new(seqid, "stop_codon", 48411, 48413,
                                             GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene2, (GtFeatureNode *)mrna2);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)exon4);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)exon5);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)exon6);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)start2);
  gt_feature_node_add_child((GtFeatureNode *)mrna2, (GtFeatureNode *)stop2);
  gt_array_add(genes, gene2);

  GtGenomeNode *gene3  = gt_feature_node_new(seqid, "gene", 88551, 92176,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *mrna3  = gt_feature_node_new(seqid, "mRNA", 88551, 92176,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon7  = gt_feature_node_new(seqid, "exon", 88551, 89029,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon8  = gt_feature_node_new(seqid, "exon", 89265, 89549,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon9  = gt_feature_node_new(seqid, "exon", 90074, 90413,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon10 = gt_feature_node_new(seqid, "exon", 90728, 90833,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon11 = gt_feature_node_new(seqid, "exon", 91150, 91362,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *exon12 = gt_feature_node_new(seqid, "exon", 91810, 92176,
                                             GT_STRAND_REVERSE);
  GtGenomeNode *start3 = gt_feature_node_new(seqid, "start_codon", 91961, 91963,
                                             GT_STRAND_FORWARD);
  GtGenomeNode *stop3  = gt_feature_node_new(seqid, "stop_codon", 88892, 88894,
                                             GT_STRAND_FORWARD);
  gt_feature_node_add_child((GtFeatureNode *)gene3, (GtFeatureNode *)mrna3);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon7);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon8);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon9);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon10);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon11);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)exon12);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)start3);
  gt_feature_node_add_child((GtFeatureNode *)mrna3, (GtFeatureNode *)stop3);
  gt_array_add(genes, gene3);

  gt_str_delete(seqid);
  return genes;
}

GtFeatureNode *agn_test_data_eden()
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
  gt_genome_node_ref(feature);
  feature = gt_feature_node_new(seqid, "exon", 1300, 1500, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna3fn, featurefn);
  feature = gt_feature_node_new(seqid, "exon", 3000, 3902, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);
  gt_genome_node_ref(feature);
  feature = gt_feature_node_new(seqid, "exon", 5000, 5500, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna2fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);
  gt_genome_node_ref(feature);
  gt_genome_node_ref(feature);
  feature = gt_feature_node_new(seqid, "exon", 7000, 9000, GT_STRAND_FORWARD);
  featurefn = (GtFeatureNode *)feature;
  gt_feature_node_add_child(mrna1fn, featurefn);
  gt_feature_node_add_child(mrna2fn, featurefn);
  gt_feature_node_add_child(mrna3fn, featurefn);
  gt_genome_node_ref(feature);
  gt_genome_node_ref(feature);

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

  gt_str_delete(seqid);
  return genefn;
}
