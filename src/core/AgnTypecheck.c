#include "extended/feature_node_iterator_api.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

bool agn_typecheck_cds(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "CDS") ||
         gt_feature_node_has_type(fn, "coding sequence") ||
         gt_feature_node_has_type(fn, "coding_sequence");
}

GtUword agn_typecheck_count(GtFeatureNode *fn, bool (*func)(GtFeatureNode *))
{
  GtUword count = 0;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *feature;
  for(feature  = gt_feature_node_iterator_next(iter);
      feature != NULL;
      feature  = gt_feature_node_iterator_next(iter))
  {
    if(func(feature))
      count++;
  }
  gt_feature_node_iterator_delete(iter);

  return count;
}

bool agn_typecheck_exon(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "exon");
}

bool agn_typecheck_gene(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "gene");
}

bool agn_typecheck_intron(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "intron");
}

bool agn_typecheck_mrna(GtFeatureNode *fn)
{ 
  return gt_feature_node_has_type(fn, "mRNA") ||
         gt_feature_node_has_type(fn, "messenger RNA") ||
         gt_feature_node_has_type(fn, "messenger_RNA");
}

bool agn_typecheck_pseudogene(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "pseudogene");
}

GtArray *agn_typecheck_select(GtFeatureNode *fn, bool (*func)(GtFeatureNode *))
{
  GtArray *children = gt_array_new( sizeof(GtFeatureNode *) );
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    if(func(current))
      gt_array_add(children, current);
  }
  gt_feature_node_iterator_delete(iter);
  gt_array_sort(children, (GtCompare)agn_genome_node_compare);
  return children;
}

bool agn_typecheck_start_codon(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "start_codon") ||
         gt_feature_node_has_type(fn, "start codon") ||
         gt_feature_node_has_type(fn, "initiation codon");
}

bool agn_typecheck_stop_codon(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "stop_codon") ||
         gt_feature_node_has_type(fn, "stop codon");
         // What about 'termination codon'?
}

bool agn_typecheck_transcript(GtFeatureNode *fn)
{
  return agn_typecheck_mrna(fn) ||
         gt_feature_node_has_type(fn, "tRNA") ||
         gt_feature_node_has_type(fn, "transfer RNA") ||
         gt_feature_node_has_type(fn, "rRNA") ||
         gt_feature_node_has_type(fn, "ribosomal RNA");
}

bool agn_typecheck_utr(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "UTR") ||
         gt_feature_node_has_type(fn, "untranslated region") ||
         gt_feature_node_has_type(fn, "untranslated_region") ||
         agn_typecheck_utr3p(fn) ||
         agn_typecheck_utr5p(fn);
}

bool agn_typecheck_utr3p(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "3' UTR") ||
         gt_feature_node_has_type(fn, "3'UTR") ||
         gt_feature_node_has_type(fn, "three prime UTR") ||
         gt_feature_node_has_type(fn, "three_prime_UTR") ||
         gt_feature_node_has_type(fn, "three prime untranslated region") ||
         gt_feature_node_has_type(fn, "three_prime_untranslated_region");
}

bool agn_typecheck_utr5p(GtFeatureNode *fn)
{
  return gt_feature_node_has_type(fn, "5' UTR") ||
         gt_feature_node_has_type(fn, "5'UTR") ||
         gt_feature_node_has_type(fn, "five prime UTR") ||
         gt_feature_node_has_type(fn, "five_prime_UTR") ||
         gt_feature_node_has_type(fn, "five prime untranslated region") ||
         gt_feature_node_has_type(fn, "five_prime_untranslated_region");
}