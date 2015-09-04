/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include "AgnLocusMapVisitor.h"
#include "AgnTypecheck.h"
#include "AgnUtils.h"

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
  bool useacc;
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
agn_locus_map_stream_new(GtNodeStream *in, FILE *genefh, FILE *mrnafh,
                         bool useacc)
{
  GtNodeVisitor *nv = agn_locus_map_visitor_new(genefh, mrnafh);
  if(useacc)
  {
    AgnLocusMapVisitor *mv = locus_map_visitor_cast(nv);
    agn_locus_map_visitor_use_accession(mv);
  }
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_locus_map_visitor_new(FILE *genefh, FILE *mrnafh)
{
  GtNodeVisitor *nv = gt_node_visitor_create(locus_map_visitor_class());
  AgnLocusMapVisitor *v = locus_map_visitor_cast(nv);
  v->genefh = genefh;
  v->mrnafh = mrnafh;
  v->useacc = false;
  return nv;
}

void agn_locus_map_visitor_use_accession(AgnLocusMapVisitor *mv)
{
  agn_assert(mv);
  mv->useacc = true;
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
  agn_assert(gt_feature_node_has_type(fn, "locus"));
  const char *locusid = gt_feature_node_get_attribute(fn, "ID");

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_typecheck_gene(current) && v->genefh != NULL)
    {
      const char *geneid = gt_feature_node_get_attribute(current, "accession");
      if(geneid == NULL || v->useacc == false)
        geneid = gt_feature_node_get_attribute(current, "ID");
      fprintf(v->genefh, "%s\t%s\n", geneid, locusid);
    }

    if(agn_typecheck_mrna(current) && v->mrnafh != NULL)
    {
      const char *mrnaid = gt_feature_node_get_attribute(current, "accession");
      if(mrnaid == NULL || v->useacc == false)
        mrnaid = gt_feature_node_get_attribute(current, "ID");
      fprintf(v->mrnafh, "%s\t%s\n", mrnaid, locusid);
    }
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}
