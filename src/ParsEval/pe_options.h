#ifndef PARSEVAL_OPTIONS
#define PARSEVAL_OPTIONS

#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

/**
 * @type FIXME
 */
enum PeOutFormat
{
  TEXTMODE,
  HTMLMODE,
  CSVMODE
};
typedef enum PeOutFormat PeOutFormat;

/**
 * @type Simple data structure for program options.
 */
struct ParsEvalOptions
{
  bool debug;
  FILE *outfile;
  char *outfilename;
  bool gff3;
  bool summary_only;
  bool graphics;
  AgnLocusPngMetadata pngdata;
  const char *refrfile;
  const char *predfile;
  const char *refrlabel;
  const char *predlabel;
  PeOutFormat outfmt;
  bool overwrite;
  const char *data_path;
  bool makefilter;
  GtArray *filters;
  bool verbose;
  int max_transcripts;
};
typedef struct ParsEvalOptions ParsEvalOptions;

void pe_free_option_memory(ParsEvalOptions *options);
int pe_parse_options(int argc, char **argv, ParsEvalOptions *options,
                     GtError *error);
void pe_print_usage(FILE *outstream);
void pe_set_option_defaults(ParsEvalOptions *options);

#endif
