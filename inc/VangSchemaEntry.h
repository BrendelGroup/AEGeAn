#ifndef VANG_SCHEMA_ENTRY
#define VANG_SCHEMA_ENTRY

#include "genometools.h"

//----------------------------------------------------------------------------//
// Type definitions
//----------------------------------------------------------------------------//
typedef struct VangSchemaEntry VangSchemaEntry;
typedef struct VangRelationExclusion VangRelationExclusion;


//----------------------------------------------------------------------------//
// Method prototypes
//----------------------------------------------------------------------------//

/**
 * Parse an entry from the given schema file
 *
 * @param[in] schemafile    pointer to a schema file
 * @returns                 an object representing the next entry in the file,
 *                          or NULL if the end of the file has been reached
 */
VangSchemaEntry *vang_schema_entry_next(FILE *schemafile);

/**
 * Free memory previously occupied by an entry object
 *
 * @param[in] entry    a schema entry
 */
void vang_schema_entry_delete(VangSchemaEntry *entry);

/**
 * Free memory previously occupied by a relation exclusion object
 *
 * @param[in] exclusion    a relation exclusion
 */
void vang_relation_exclusion_delete(VangRelationExclusion *exclusion);

#endif
