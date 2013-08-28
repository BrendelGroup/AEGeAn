#ifndef PE_OPTIONS
#define PE_OPTIONS

#include "genometools.h"
#include "AgnComparEval.h"

/**
 * This struct contains ParsEval's command-line options.
 */
typedef struct
{
  bool debug;
  FILE *outfile;
  char *outfilename;
  bool verbose;
  bool gff3;
  int complimit;
  bool summary_only;
  bool vectors;
  bool locus_graphics;
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
  const char *outfmt;
  bool overwrite;
  bool html;
  const char *data_path;
  bool usefilter;
  const char *filterfile;
  AgnCompareFilters filters;
  int numprocs;
  int trans_per_locus;
} PeOptions;

/**
 * Parse command-line options from the provided arguments.
 *
 * @param[in] argc    the number of arguments
 * @param[in] argv    the list of command-line arguments
 */
int pe_parse_options(int argc, char * const argv[], PeOptions *options);

/**
 * Print usage statement
 */
void pe_print_usage();

/**
 * Initialize command-line options to default values.
 */
void pe_set_option_defaults(PeOptions *options);

#endif
