#ifndef AGN_UTILS
#define AGN_UTILS

#include "core/array_api.h"
#include "extended/genome_node_api.h"

/**
 * @module AgnUtils
 *
 * Collection of assorted functions that are otherwise unrelated.
 */ //;

/**
 * @function FIXME
 */
double agn_calc_splice_complexity(GtArray *transcripts);

/**
 * @function Compare function for data type ``GtGenomeNode **``, needed for
 * sorting ``GtGenomeNode *`` stored in ``GtArray`` objects.
 */
int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b);

#endif
