#ifndef AEGEAN_TEST_DATA
#define AEGEAN_TEST_DATA

/**
 * @module AgnTestData
 * A collection of functions facilitating unit testing of various AEGeAn classes
 * and modules.
 */ //;
 
#include "genometools.h"

/**
 * @function Example from grape..
 */
GtArray *agn_test_data_grape();

/**
 * @function Example from grape: gene structure annotated with exon and start /
 * stop codon features--CDS is implicitly defined by these features.
 */
GtArray *agn_test_data_grape_codons();

/**
 * @function Create the canonical gene structure (from the GFF3 specification)
 * in memory.
 */
GtFeatureNode *agn_test_data_eden();

#endif
