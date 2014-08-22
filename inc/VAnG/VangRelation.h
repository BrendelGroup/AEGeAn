/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef VANG_RELATION
#define VANG_RELATION

#include "VangDegreeConstraint.h"

//----------------------------------------------------------------------------//
// Enumerated types and type definitions
//----------------------------------------------------------------------------//
enum SpatialConstraintEnum {SC_NONE,SC_CONTAINS,SC_WITHIN,SC_OVERLAP,SC_EXACT};
typedef struct VangRelation VangRelation;
typedef enum SpatialConstraintEnum SpatialConstraint;


//----------------------------------------------------------------------------//
// Method prototypes
//----------------------------------------------------------------------------//

/**
 * Allocate memory for a relation object
 *
 * @param[in] id                    unique identifier for this relation
 * @param[in] nodetype              the data type of the external feature to
                                    which this relation refers
 * @returns                         a relation object which can be used for
                                    feature validation
 */
VangRelation *vang_relation_new(const char *id, const char *nodetype);

/**
 * Free memory previous occupied by a relation object
 *
 * @param[in] rel    a relation object
 */
void vang_relation_delete(VangRelation *rel);

/**
 * Specify the degree constraint of the relation; for example, the relation
 *
 *    Relation=r1;Nodetype=CDS;Degree=in e 1
 *
 * declared for data type 'mRNA' dictates that the number of edges pointing to
 * any given mRNA feature from CDS features must be exactly 1, or else
 * validation fails.
 *
 * @param[out] rel           the relation object to which this degree
 *                           constraint corresponds
 * @param[in]  context       specifies whether this constraint is applied to
 *                           the indegree or outdegree
 * @param[in]  operator      specifies the operator that will be used to test /
 *                           enforce this constraint
 * @param[in]  constraint    numeric constraint on the degree
 */
void vang_relation_set_degree(VangRelation *rel, DegreeContext context,
                              DegreeOperator operator, unsigned int degree);

/**
 * Specify the attribute key used to declare this relation; default is 'Parent'
 *
 * @param[out] rel    the relation object
 * @param[in]  key    the key used to declare this relation
 */
void vang_relation_set_key(VangRelation *rel, const char *key);

/**
 * Spatial constraints can be applied to relations, although the default is to
 * apply no constraint
 *
 * @param[out] rel        the relation object
 * @param[in]  spatial    the spatial constraint to be applied to this
 *                        relation
 */
void vang_relation_set_spatial(VangRelation *rel, SpatialConstraint spatial);

/**
 * Attach a descriptive free-text note to this relation
 *
 * @param[out] rel     the relation object
 * @param[in]  note    the note text
 */
void vang_relation_set_note(VangRelation *rel, const char *note);

/**
 * Print a string representation of this relation
 *
 * @param[in]  rel          relation object
 * @param[out] outstream    output stream to which data will be printed
 */
void vang_relation_to_string(VangRelation *rel, FILE *outstream);

/**
 *
 */
SpatialConstraint vang_spatial_constraint_parse(const char *spatial_string);

#endif

