/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_GENE_LOCUS_MAPPING
#define AEGEAN_GENE_LOCUS_MAPPING

#include "core/error_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"
#include "AgnLocus.h"

/**
 * @class AgnGeneLocusMapping
 *
 * Class for maintaining a mapping of gene IDs to locus IDs, working with the
 * mapping in memory, and writing the mapping to disk.
 */
typedef struct AgnGeneLocusMapping AgnGeneLocusMapping;

/**
 * @function Inspect a locus object and add all of its associated genes to the
 * map.
 */
void agn_gene_locus_mapping_add(AgnGeneLocusMapping *map, AgnLocus *locus);

/**
 * @function Write the map to a file and call the class destructor to free
 * memory.
 */
void agn_gene_locus_mapping_close(AgnGeneLocusMapping *map);

/**
 * @function Destructor.
 */
void agn_gene_locus_mapping_delete(AgnGeneLocusMapping *map);

/**
 * @function Return a list of all the gene IDs associated with the given locus
 * position.
 */
GtStrArray *
agn_gene_locus_mapping_get_geneids_for_locus(AgnGeneLocusMapping *map,
                                             const char *locuspos,
                                             GtError *error);

/**
 * @function Given a gene ID, return the associated locus ID.
 */
GtStr *agn_gene_locus_mapping_get_locus(AgnGeneLocusMapping *map,
                                        const char *geneid);

/**
 * @function Read a gene-locus mapping from a file.
 */
AgnGeneLocusMapping *agn_gene_locus_mapping_open(const char *filepath);

/**
 * @function Create a new, empty gene-locus mapping in memory.
 */
AgnGeneLocusMapping *agn_gene_locus_mapping_new(const char *filepath);

/**
 * @function Remove a gene ID from the map and return the associated locus ID.
 * Caller is responsible to free the locus ID.
 */
GtStr *agn_gene_locus_mapping_unmap_gene(AgnGeneLocusMapping *map,
                                         const char *geneid);

#endif
