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
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *locus_map_visitor_class();

/**
 * @function FIXME
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
  return gt_visitor_stream_new(in, nv);
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
  agn_assert(gt_feature_node_has_type(fn, "locus"));
  const char *locuslabel = agn_feature_node_get_label(fn);

  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
  GtFeatureNode *current;
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_typecheck_gene(current) && v->genefh != NULL)
    {
      const char *genelabel = agn_feature_node_get_label(current);
      fprintf(v->genefh, "%s\t%s\n", genelabel, locuslabel);
    }

    if(agn_typecheck_mrna(current) && v->mrnafh != NULL)
    {
      const char *mrnalabel = agn_feature_node_get_label(current);
      fprintf(v->mrnafh, "%s\t%s\n", mrnalabel, locuslabel);
    }
  }
  gt_feature_node_iterator_delete(iter);
  return 0;
}
