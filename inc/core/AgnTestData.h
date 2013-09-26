#ifndef AEGEAN_TEST_DATA
#define AEGEAN_TEST_DATA

/**
 * @module AgnTestData
 * A collection of functions facilitating unit testing of various AEGeAn classes
 * and modules.
 */ //;
 
 #include "genometools.h"

/**
 * @function Create an array containing 3 gene features to be used for unit
 * testing.
 */
GtArray *agn_test_data_genes_codons();

/**
 * @function Create the canonical gene structure (from the GFF3 specification)
 * in memory.
 */
GtFeatureNode *agn_test_data_eden();

#endif
