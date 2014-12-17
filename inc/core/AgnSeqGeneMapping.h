/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#ifndef AEGEAN_SEQ_GENE_MAPPING
#define AEGEAN_SEQ_GENE_MAPPING

#include "core/str_array_api.h"
#include "AgnLocus.h"

/**
 * @class AgnSeqGeneMapping
 *
 * Class for maintaining a mapping of sequence IDs to lists of corresponding
 * gene IDs, working with the mapping in memory, and writing the mapping to
 * disk.
 */
typedef struct AgnSeqGeneMapping AgnSeqGeneMapping;

/**
 * @function Inspect a locus object and add all of its associated genes to the
 * map.
 */
void agn_seq_gene_mapping_add(AgnSeqGeneMapping *map, AgnLocus *locus);

/**
 * @function Write the map to a file and call the class destructor to free
 * memory.
 */
void agn_seq_gene_mapping_close(AgnSeqGeneMapping *map);

/**
 * @function Destructor.
 */
void agn_seq_gene_mapping_delete(AgnSeqGeneMapping *map);

/**
 * @function Return an array of gene IDs associated with the specified sequence,
 * or NULL if there are none.
 */
GtStrArray *agn_seq_gene_mapping_get_genes(AgnSeqGeneMapping *map,
                                           const char *seqid);

/**
 * @function Test whether a sequence is present in the map.
 */
bool agn_seq_gene_mapping_has_seqid(AgnSeqGeneMapping *map, const char *seqid);

/**
 * @function Read a sequence-gene mapping from a file.
 */
AgnSeqGeneMapping *agn_seq_gene_mapping_open(const char *filepath);

/**
 * @function Create a new, empty sequence-gene mapping in memory.
 */
AgnSeqGeneMapping *agn_seq_gene_mapping_new(const char *filepath);

/**
 * @function Search the genes associated with the given sequence, and if
 * `geneid` is found, remove it and return true; otherwise return false.
 */
bool agn_seq_gene_mapping_unmap_gene(AgnSeqGeneMapping *map,
                                     const char *seqid,
                                     const char *geneid);

/**
 * @function Remove the specified sequence from the map and return an array of
 * any associated gene IDs.
 */
GtStrArray *agn_seq_gene_mapping_unmap_seqid(AgnSeqGeneMapping *map,
                                             const char *seqid);

#endif
