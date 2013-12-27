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
 * @function Returns true if the given feature is an exon; false otherwise.
 */
bool agn_typecheck_exon(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is an intron; false otherwise.
 */
bool agn_typecheck_intron(GtFeatureNode *fn);

/**
 * @function Returns true if the given feature is an mRNA; false otherwise.
 */
bool agn_typecheck_mrna(GtFeatureNode *fn);

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
