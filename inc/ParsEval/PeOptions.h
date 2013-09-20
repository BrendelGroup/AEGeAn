#ifndef PE_OPTIONS
#define PE_OPTIONS

/**
 * @module PeOptions
 * Module for defining ParsEval's command line options, parsing them from the
 * command line, and storing/accessing them during runtime.
 */ //;

#include "genometools.h"
#include "AgnComparEval.h"

/**
 * @type This struct defines ParsEval's command-line options.
 */
struct PeOptions
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
  const char *data_path;
  bool makefilter;
  bool usefilter;
  const char *filterfile;
  AgnCompareFilters filters;
  int numprocs;
  int trans_per_locus;
};
typedef struct PeOptions PeOptions;

/**
 * @function Parse command-line options from the provided arguments.
 *
 * @param[in] argc    the number of arguments
 * @param[in] argv    the list of command-line arguments
 */
int pe_parse_options(int argc, char * const argv[], PeOptions *options);

/**
 * @function Print usage statement.
 */
void pe_print_usage();

/**
 * @function Initialize command-line options to default values.
 */
void pe_set_option_defaults(PeOptions *options);

#endif
