#ifndef PE_REPORTS
#define PE_REPORTS

/**
 * @module PeReports
 * Module for generating ParsEval's reports in text, CSV, or HTML format.
 */

#include "genometools.h"
#include "AgnGeneLocus.h"
#include "AgnLocusIndex.h"
#include "AgnLogger.h"
#include "PeOptions.h"

#define PE_GENE_LOCUS_GRAPHIC_MIN_WIDTH 650

/**
 * @function Get the filename for printing this locus' results.
 */
void pe_gene_locus_get_filename(AgnGeneLocus *locus, char *buffer,
                                const char *dirpath);

/**
 * @function Determine the optimal width (in pixels) of the PNG graphic used to
 * visualize this locus.
 */
GtUword pe_gene_locus_get_graphic_width(AgnGeneLocus *locus);

/**
 * @function Get the filename for printing this locus' graphic, write it to
 * ``buffer``.
 */
void pe_gene_locus_get_png_filename(AgnGeneLocus *locus, char *buffer,
                                    const char *dirpath);

/**
 * @function Generate a string containing the current time. User is responsible
 * to free the string.
 */
char *pe_get_start_time();

/**
 * @function Print the first row of the CSV output file, containing column
 * headers.
 */
void pe_print_csv_header(FILE *outstream);

/**
 * @function Write the report for this locus.
 */
void pe_gene_locus_print_results(AgnGeneLocus *locus, FILE *outstream,
                                 PeOptions *options);

/**
 * @function Write the report for this locus in CSV format.
 */
void pe_gene_locus_print_results_csv(AgnGeneLocus *locus, FILE *outstream,
                                     PeOptions *options);

/**
 * @function Write the report for this locus in HTML format.
 */
void pe_gene_locus_print_results_html(AgnGeneLocus *locus, PeOptions *options);

/**
 * @function Print footer for HTML files.
 */
void pe_print_html_footer(FILE *outstream);

/**
 * @function In the HTML output, there is a list of loci for each sequence. This
 * function prints a line in that file.
 */
void pe_print_locus_to_seqfile(FILE *seqfile, GtUword start, GtUword end,
                               GtUword length, GtUword refr_transcripts,
                               GtUword pred_transcripts,
                               AgnCompSummary *comparisons);

/**
 * @function Print header for sequence-specific HTML files.
 */
void pe_print_seqfile_header(FILE *outstream, const char *seqid);

/**
 * @function Print footer for sequence-specific HTML files.
 */
void pe_print_seqfile_footer(FILE *outstream);

/**
 * @function Print the ParsEval summary. FIXME this could use refactoring
 */
void pe_print_summary(const char *start_time, int argc, char * const argv[],
                      GtStrArray *seqids, AgnCompEvaluation *summary_data,
                      GtArray *seq_summary_data, FILE *outstream,
                      PeOptions *options);
void pe_print_summary_html(const char *start_time, int argc,
                           char * const argv[], GtStrArray *seqids,
                           AgnCompEvaluation *summary_data,
                           GtArray *seq_summary_data, FILE *outstream,
                           PeOptions *options);

/**
 * @function ParsEval uses sequence IDs to create temporary and permanent output
 * files. Sequence IDs that contain characters that are not supported for POSIX
 * filenames cause problems. This function checks for those problems.
 *
 * @param[in]  seqid     the sequence ID to be checked
 * @param[out] logger    object to which an error message will be written if the
 *                       seqid is invalid
 */
void pe_seqid_check(const char *seqid, AgnLogger *logger);

/**
 * @function Comparison function used to set the track order in PNG graphics.
 */
int pe_track_order(const char *s1, const char *s2, void *data);

#endif
