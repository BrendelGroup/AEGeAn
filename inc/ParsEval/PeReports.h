#ifndef PE_REPORTS
#define PE_REPORTS

#include "genometools.h"
#include "AgnPairwiseCompareLocus.h"
#include "PeComparEval.h"
#include "PeOptions.h"

#define AGN_PAIRWISE_COMPARE_LOCUS_GRAPHIC_MIN_WIDTH 650

/**
* This data structure provides a convenient container for metadata needed to
* produce a locus PNG graphic.
*/
typedef struct
{
  char filename[512];
  char stylefile[512];
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
  unsigned long graphic_width;
  int (*track_order_func)(const char *s1, const char *s2, void *data);
} AgnPairwiseCompareLocusPngMetadata;

/**
* This data structure provides a summary of the data and comparisons associated
* with a given locus.
*/
typedef struct
{
  unsigned long start;
  unsigned long end;
  unsigned long length;
  unsigned long refr_transcripts;
  unsigned long pred_transcripts;
  unsigned long reported;
  unsigned long total;
  AgnCompSummary counts;
} AgnPairwiseCompareLocusSummary;

/**
 * Add the locus' comparison statistics to a set of aggregate statistics.
 *
 * @param[in]  locus    the locus annotation
 * @param[out] data     summary counts, stats, and results to which the locus
 *                      data will be aggregated
 */
void agn_pairwise_compare_locus_aggregate_results
(
  AgnPairwiseCompareLocus *locus,
  PeCompEvaluation *data
);

#ifndef WITHOUT_CAIRO
/**
 * Print a PNG graphic for this locus.
 *
 * @param[in] locus    the locus
 */
void agn_pairwise_compare_locus_print_png
(
  AgnPairwiseCompareLocus *locus,
  AgnPairwiseCompareLocusPngMetadata *metadata
);
#endif

/**
 * Get the filename for printing this locus' results.
 *
 * @param[in]  locus     the locus
 * @param[out] buffer    the string to which the filename will be written
 */
void pe_gene_locus_get_filename( AgnPairwiseCompareLocus *locus, char *buffer,
                                 const char *dirpath );

/**
 * Determine the width of the PNG graphic used to visualize this locus
 *
 * @param[in]    the locus
 * @returns      the width (in pixels) to use
 */
unsigned long pe_gene_locus_get_graphic_width(AgnPairwiseCompareLocus *locus);

/**
 * Get the filename for printing this locus' graphic.
 *
 * @param[in]  locus       the locus
 * @param[out] filename    string to which filename will be written
 */
void pe_gene_locus_get_png_filename( AgnPairwiseCompareLocus *locus,
                                     char *buffer, const char *dirpath );

/**
 * Write the report for this locus.
 *
 * @param[in]  locus        the locus
 * @param[out] outstream    the file to which the report will be written
 */
void pe_gene_locus_print_results( AgnPairwiseCompareLocus *locus,
                                  FILE *outstream, PeOptions *options );

/**
 * Write the report for this locus in CSV format.
 *
 * @param[in]  locus        the locus
 * @param[out] outstream    the file to which the report will be written
 */
void pe_gene_locus_print_results_csv( AgnPairwiseCompareLocus *locus,
                                      FILE *outstream, PeOptions *options );

/**
 * Write the report for this locus in HTML format.
 *
 * @param[in] locus    the locus
 */
void pe_gene_locus_print_results_html( AgnPairwiseCompareLocus *locus,
                                       PeOptions *options );

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
void pe_print_locus_to_seqfile( FILE *seqfile, unsigned long start,
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
 * Print the ParsEval summary.
 */
void pe_print_summary( const char *start_time, int argc, char * const argv[],
                       GtStrArray *seqids, PeCompEvaluation *summary_data,
                       PeCompEvaluation *seq_summary_data, FILE *outstream,
                       PeOptions *options );
void pe_print_summary_html( const char *start_time, int argc,
                            char * const argv[], GtStrArray *seqids,
                            PeCompEvaluation *summary_data,
                            PeCompEvaluation *seq_summary_data, FILE *outstream,
                            PeOptions *options );

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
