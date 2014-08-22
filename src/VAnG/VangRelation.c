/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <string.h>
#include "genometools.h"
#include "VangRelation.h"

//----------------------------------------------------------------------------//
// Private data structure definition
//----------------------------------------------------------------------------//
struct VangRelation
{
  char *id;
  DegreeConstraint dc;
  char *nodetype;
  char *key;
  SpatialConstraint spatial;
  char *note;
};


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//
VangRelation *vang_relation_new(const char *id, const char *nodetype)
{
  VangRelation *rel = gt_malloc( sizeof(VangRelation) );
  rel->id = gt_cstr_dup(id);
  rel->nodetype = gt_cstr_dup(nodetype);

  // Apply default degree constraint: indegree == 1; see VangRelation.h for a
  // more detailed explanation
  rel->dc.context = DEGREE_IN;
  rel->dc.operator = EQUALS;
  rel->dc.degree = 1;

  // Apply default key (Parent), spatial constraint (none), and note (undefined)
  rel->key = gt_cstr_dup("Parent");
  rel->spatial = SC_NONE;
  rel->note = NULL;

  return rel;
};

void vang_relation_delete(VangRelation *rel)
{
  gt_free(rel->id);
  gt_free(rel->nodetype);
  gt_free(rel->key);
  if(rel->note != NULL)
    gt_free(rel->note);
  gt_free(rel);
  rel = NULL;
}

void vang_relation_set_degree(VangRelation *rel, DegreeContext context,
                              DegreeOperator operator, unsigned int degree)
{
  rel->dc.context = context;
  rel->dc.operator = operator;
  rel->dc.degree = degree;
}

void vang_relation_set_key(VangRelation *rel, const char *key)
{
  gt_free(rel->key);
  rel->key = gt_cstr_dup(key);
}

void vang_relation_set_spatial(VangRelation *rel, SpatialConstraint spatial)
{
  rel->spatial = spatial;
}

void vang_relation_set_note(VangRelation *rel, const char *note)
{
  if(rel->note != NULL)
    gt_free(rel->note);
  rel->note = gt_cstr_dup(note);
}

void vang_relation_to_string(VangRelation *rel, FILE *outstream)
{
  const char *spatial;
  switch(rel->spatial)
  {
    case SC_NONE:
      spatial = "none";
      break;
    case SC_CONTAINS:
      spatial = "contains";
      break;
    case SC_WITHIN:
      spatial = "within";
      break;
    case SC_OVERLAP:
      spatial = "overlap";
      break;
    case SC_EXACT:
      spatial = "exact";
      break;
    default:
      spatial = "undef";
      break;
  }

  fprintf(outstream, "Relation=%s;Nodetype=%s;Degree=", rel->id, rel->nodetype);
  vang_degree_constraint_to_string(&rel->dc, outstream);
  fprintf(outstream, ";Key=%s;Spatial=%s", rel->key, spatial);
  if(rel->note != NULL)
    fprintf(outstream, ";Note=%s", rel->note);
}

SpatialConstraint vang_spatial_constraint_parse(const char *spatial_string)
{
  if(strcmp(spatial_string, "none") == 0)
    return SC_NONE;
  else if(strcmp(spatial_string, "contains") == 0)
    return SC_CONTAINS;
  else if(strcmp(spatial_string, "within") == 0)
    return SC_WITHIN;
  else if(strcmp(spatial_string, "overlap") == 0)
    return SC_OVERLAP;
  else if(strcmp(spatial_string, "exact") == 0)
    return SC_EXACT;
  else
  {
    fprintf(stderr, "error: unknown spatial constraint '%s'\n", spatial_string);
    exit(1);
  }

  return SC_NONE;
}

