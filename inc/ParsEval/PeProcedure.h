#ifndef PE_PROCEDURE
#define PE_PROCEDURE

/**
 * @module PeProcedure
 * This module implements the main steps of the ParsEval program.
 */

#include <string.h>
#include "genometools.h"
#include "AgnComparEval.h"
#include "AgnLocusIndex.h"
#include "AgnLogger.h"
#include "AgnUtils.h"
#include "PeOptions.h"

/**
 * @type Supplemental data structure used during comparative analysis of gene
 * loci.
 */
struct PeAnalysisData
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;
  PeOptions *options;
  FILE *seqfile;
};
typedef struct PeAnalysisData PeAnalysisData;

/**
 * @function Aggregate locus-specific comparison statistics at the sequence level
 * and over all sequences.
 */
void pe_aggregate_results(AgnCompEvaluation *overall_eval,
                          GtArray **seqlevel_evalsp, GtArray *loci,
                          GtArray *seqfiles, GtHashmap *comp_evals,
                          GtHashmap *locus_summaries, PeOptions *options);

/**
 * @function Perform comparison of annotations (reference vs prediction) for
 * each locus in the given locus index.
 */
void pe_comparative_analysis(AgnLocusIndex *locusindex, GtHashmap **comp_evalsp,
                             GtHashmap **locus_summariesp, GtStrArray *seqids,
                             GtArray *seqfiles, GtArray *loci,
                             PeOptions *options);

/**
 * @function Load gene annotations into memory and identify gene loci.
 */
GtUword pe_load_and_parse_loci(AgnLocusIndex **locusindexp, GtArray **locip,
                               GtStrArray **seqidsp, PeOptions *options,
                               AgnLogger *logger);

/**
 * @function Collect information from the given locus following comparative
 * analysis.
 */
void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

/**
 * @function Collect information from the given locus prior to comparative
 * analysis.
 */
void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data);

/**
 * @function Prepare output files.
 */
GtArray *pe_prep_output(GtStrArray *seqids, PeOptions *options);

/**
 * @function Print the ParsEval summary and, if necessary, combine temporary
 * output files to create the final output.
 */
void pe_print_combine_output(GtStrArray *seqids, GtArray *seqfiles,
                             PeOptions *options);

#endif
