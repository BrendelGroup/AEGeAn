/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
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
 * Get the data type associated with this entry
 *
 * @param[in] entry    a schema entry
 * @returns            the entry's data type
 */
const char *vang_schema_entry_get_datatype(VangSchemaEntry *entry);

/**
 * Print a string representation of this schema entry
 *
 * @param[in]  entry        schema entry object
 * @param[out] outstream    output stream to which data will be printed
 */
void vang_schema_entry_to_string(VangSchemaEntry *entry, FILE *outstream);

/**
 * Free memory previously occupied by a relation exclusion object
 *
 * @param[in] exclusion    a relation exclusion
 */
void vang_relation_exclusion_delete(VangRelationExclusion *exclusion);

#endif
