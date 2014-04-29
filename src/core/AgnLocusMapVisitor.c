#include "AgnLocusMapVisitor.h"
#include "AgnTypecheck.h"

#define locus_map_visitor_cast(GV)\
        gt_node_visitor_cast(locus_map_visitor_class(), GV)

//----------------------------------------------------------------------------//
// Data structure definition
//----------------------------------------------------------------------------//

struct AgnLocusMapVisitor
{
  const GtNodeVisitor parent_instance;
  FILE *genefh;
  FILE *mrnafh;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *locus_map_visitor_class();

/**
 * @function For any gene feature with attribute 'pseudo=true', set type to
 * 'pseudogene'.
 */
static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream*
agn_locus_map_stream_new(GtNodeStream *in, FILE *genefh, FILE *mrnafh)
{
  GtNodeVisitor *nv = agn_locus_map_visitor_new(genefh, mrnafh);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_locus_map_visitor_new(FILE *genefh, FILE *mrnafh)
{
  GtNodeVisitor *nv = gt_node_visitor_create(locus_map_visitor_class());
  AgnLocusMapVisitor *v = locus_map_visitor_cast(nv);
  v->genefh = genefh;
  v->mrnafh = mrnafh;
  return nv;
}

static const GtNodeVisitorClass *locus_map_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnLocusMapVisitor), NULL, NULL,
                                    visit_feature_node, NULL, NULL, NULL);
  }
  return nvc;
}

static int
visit_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn, GtError *error)
{
  AgnLocusMapVisitor *v = locus_map_visitor_cast(nv);
  gt_error_check(error);
  gt_assert(gt_feature_node_has_type(fn, "locus"));
  char lid[1024];
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)fn);
  GtRange range = gt_genome_node_get_range((GtGenomeNode *)fn);
  sprintf(lid, "%s:%lu-%lu", gt_str_get(seqid), range.start, range.end);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_typecheck_gene(current) && v->genefh != NULL)
    {
      const char *geneid = gt_feature_node_get_attribute(current, "ID");
      fprintf(v->genefh, "%s\t%s\n", geneid, lid);
    }

    if(agn_typecheck_mrna(current) && v->mrnafh != NULL)
    {
      const char *mrnaid = gt_feature_node_get_attribute(current, "ID");
      fprintf(v->mrnafh, "%s\t%s\n", mrnaid, lid);
    }
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}
