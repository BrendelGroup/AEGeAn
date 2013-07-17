#ifndef PE_REPORTS
#define PE_REPORTS

#include "genometools.h"
#include "AgnGeneLocus.h"
#include "AgnLocusIndex.h"
#include "AgnLogger.h"
#include "PeComparEval.h"
#include "PeOptions.h"

#define PE_GENE_LOCUS_GRAPHIC_MIN_WIDTH 650

/**
 * FIXME
 */
void pe_agg_results(PeCompEvaluation *overall_eval, GtArray **seqlevel_evalsp,
                    GtArray *loci, GtArray *seqfiles, GtHashmap *comp_evals,
                    GtHashmap *locus_summaries, PeOptions *options);

/**
 * Add the locus' comparison statistics to a set of aggregate statistics.
 *
 * @param[in]  locus    the locus annotation
 * @param[out] data     summary counts, stats, and results to which the locus
 *                      data will be aggregated
 */
void pe_gene_locus_aggregate_results(AgnGeneLocus *locus,
                                     PeCompEvaluation *data);

/**
 * Get the filename for printing this locus' results.
 *
 * @param[in]  locus     the locus
 * @param[out] buffer    the string to which the filename will be written
 */
void pe_gene_locus_get_filename(AgnGeneLocus *locus, char *buffer,
                                const char *dirpath );

/**
 * Determine the width of the PNG graphic used to visualize this locus
 *
 * @param[in]    the locus
 * @returns      the width (in pixels) to use
 */
unsigned long pe_gene_locus_get_graphic_width(AgnGeneLocus *locus);

/**
 * Get the filename for printing this locus' graphic.
 *
 * @param[in]  locus       the locus
 * @param[out] filename    string to which filename will be written
 */
void pe_gene_locus_get_png_filename(AgnGeneLocus *locus, char *buffer,
                                    const char *dirpath );

/**
 * FIXME
 */
char *pe_get_start_time();

/**
 * FIXME
 */
void pe_print_csv_header(FILE *outstream);

/**
 * Write the report for this locus.
 *
 * @param[in]  locus        the locus
 * @param[out] outstream    the file to which the report will be written
 */
void pe_gene_locus_print_results(AgnGeneLocus *locus, FILE *outstream,
                                 PeOptions *options );

/**
 * Write the report for this locus in CSV format.
 *
 * @param[in]  locus        the locus
 * @param[out] outstream    the file to which the report will be written
 */
void pe_gene_locus_print_results_csv(AgnGeneLocus *locus, FILE *outstream,
                                     PeOptions *options);

/**
 * Write the report for this locus in HTML format.
 *
 * @param[in] locus    the locus
 */
void pe_gene_locus_print_results_html(AgnGeneLocus *locus, PeOptions *options);

/**
 * FIXME
 */
unsigned long pe_load_and_parse_loci(AgnLocusIndex **locusindexp,
                                     GtArray **locip, GtStrArray **seqidsp,
                                     PeOptions *options, AgnLogger *logger);

/**
 * FIXME
 */
GtArray *pe_prep_output(GtStrArray *seqids, PeOptions *options);

/**
 * FIXME
 */
void pe_print_combine_output(GtStrArray *seqids, GtArray *seqfiles,
                             PeOptions *options);

/**
 * Print footer for HTML files.
 *
 * @param[out] outstream    the output stream to which the HTML will be written
 */
void pe_print_html_footer(FILE *outstream);

/**
 * In the HTML output, there is a list of loci for each sequence. This function
 * prints a line in that file.
 *
 * @param[out] seqfile             output stream to which the HTML will be
 *                                 written
 * @param[in]  start               start coordinate of the locus
 * @param[in]  end                 end coordinate of the locus
 * @param[in]  length              the locus length
 * @param[in]  refr_transcripts    the number of reference transcripts for the
 *                                 locus
 * @param[in]  pred_transcripts    the number of prediction transcripts for the
 *                                 locus
 * @param[in]  comparisons         a breakdown of then comparisons for the locus
 *                                 by results
 */
void pe_print_locus_to_seqfile(FILE *seqfile, unsigned long start,
                               unsigned long end, unsigned long length,
                               unsigned long refr_transcripts,
                               unsigned long pred_transcripts,
                               AgnCompSummary *comparisons );

/**
 * Print header for sequence-specific HTML files.
 *
 * @param[out] outstream    the output stream to which the HTML will be written
 * @param[in]  seqid        the sequence ID
 */
void pe_print_seqfile_header(FILE *outstream, const char *seqid);

/**
 * Print footer for sequence-specific HTML files.
 *
 * @param[out] outstream    the output stream to which the HTML will be written
 */
void pe_print_seqfile_footer(FILE *outstream);

/**
 * Print the ParsEval summary. FIXME this could use refactoring
 *
 * @param[in] start_time          time at which ParsEval was launched
 * @param[in] argc                number of command-line arguments
 * @param[in] argv                vector of command-line arguments
 * @param[in] seqids              list of sequence IDs
 * @param[in] summary_data        overall comparison evaluations
 * @param[in] seq_summary_data    sequence-level comparison evaluations
 * @param[in] outstream           file stream to which summary will be written
 * @param[in] options             ParsEval parameter values
 */
void pe_print_summary(const char *start_time, int argc, char * const argv[],
                      GtStrArray *seqids, PeCompEvaluation *summary_data,
                      GtArray *seq_summary_data, FILE *outstream,
                      PeOptions *options);
void pe_print_summary_html(const char *start_time, int argc,
                           char * const argv[], GtStrArray *seqids,
                           PeCompEvaluation *summary_data,
                           GtArray *seq_summary_data, FILE *outstream,
                           PeOptions *options);

/**
 * FIXME
 */
void pe_seqid_check(const char *seqid, AgnLogger *logger);

/**
 * Comparison function used to set the track order in PNG graphics.
 *
 * @param[in] s1    string representing a track
 * @param[in] s2    string representing another track
 * @returns         1 if s1 comes before (above) s2, -1 if s2 comes before s1,
 *                  and 0 if it does not matter
 */
int pe_track_order(const char *s1, const char *s2, void *data);

#endif
