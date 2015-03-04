/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef AEGEAN_TYPECHECK
#define AEGEAN_TYPECHECK

#include "extended/feature_node_api.h"

/**
 * @module AgnTypecheck
 *
 * Functions for testing feature types.
 */ //;

/**
 * @function Returns true if the given feature is a CDS; false otherwise.
 */
bool agn_typecheck_cds(GtFeatureNode *fn);

/**
 * @function Count the number of ``fn``'s children that have the given type.
 */
GtUword agn_typecheck_count(GtFeatureNode *fn, bool (*func)(GtFeatureNode *));

/**
 * @function Returns true if the given feature is an exon; false otherwise.
 */
bool agn_typecheck_exon(GtFeatureNode *fn);

/**
 * @function Traverse the feature graph starting at `root` and add up the length
 * of all features matching the given selection function `func`.
 */
GtUword agn_typecheck_feature_combined_length(GtFeatureNode *root,
                                              bool (*func)(GtFeatureNode *));

/**
 * @function Returns true if the given feature is a gene; false otherwise.
 */
bool agn_typecheck_gene(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is an intron; false otherwise.
 */
bool agn_typecheck_intron(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is an mRNA; false otherwise.
 */
bool agn_typecheck_mrna(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is declared as a pseudogene;
 * false otherwise.
 */
bool agn_typecheck_pseudogene(GtFeatureNode *fn);

/**
 * @function Gather the children of a given feature that have a certain type.
 * Type is tested by ``func``, which accepts a single ``GtFeatureNode`` object.
 */
GtArray *agn_typecheck_select(GtFeatureNode *fn, bool (*func)(GtFeatureNode *));

/**
 * @function Returns true if the given feature is a start codon; false
 * otherwise.
 */
bool agn_typecheck_start_codon(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is a stop codon; false otherwise.
 */
bool agn_typecheck_stop_codon(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is an mRNA, tRNA, or rRNA; false
 * otherwise.
 */
bool agn_typecheck_transcript(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is a UTR; false otherwise.
 */
bool agn_typecheck_utr(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is a 3' UTR; false otherwise.
 */
bool agn_typecheck_utr3p(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is a 5' UTR; false otherwise.
 */
bool agn_typecheck_utr5p(GtFeatureNode *fn);

#endif
