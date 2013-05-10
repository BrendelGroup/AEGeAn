#include "AgnLocus.h"

void agn_locus_add(AgnLocus *locus, GtFeatureNode *gene)
{
  gt_dlist_add(locus->genes, gene);
  GtRange generange = gt_genome_node_get_range((GtGenomeNode *)gene);
  if(locus->range.start == 0 || generange.start < locus->range.start)
    locus->range.start = generange.start;
  if(locus->range.end == 0 || generange.end > locus->range.end)
    locus->range.end = generange.end;
}

void agn_locus_delete(AgnLocus *locus)
{
  gt_dlist_delete(locus->genes);
  gt_free(locus);
  locus = NULL;
}

AgnLocus *agn_locus_new(const char *seqid)
{
  AgnLocus *locus = gt_malloc( sizeof(AgnLocus) );
  locus->seqid = seqid;
  locus->genes = gt_dlist_new((GtCompare)gt_genome_node_cmp);
  locus->range.start = 0;
  locus->range.end = 0;

  return locus;
}

void agn_locus_print(AgnLocus *locus, FILE *outstream)
{
  fprintf( outstream,
           "%s\tgeneannology:locus\tlocus\t%lu\t%lu\t.\t.\t.\tnum_genes=%lu\n",
           locus->seqid, locus->range.start, locus->range.end,
           gt_dlist_size(locus->genes) );
}

void agn_locus_stringify(AgnLocus *locus, char *string)
{
  sprintf( string, "%s_%lu-%lu", locus->seqid, locus->range.start,
           locus->range.end );
}
