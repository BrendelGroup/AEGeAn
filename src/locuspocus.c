#include <omp.h>
#include <getopt.h>
#include "genometools.h"
#include "AgnLocus.h"
#include "AgnLocusIndex.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  int numprocs;
  FILE *outstream;
  bool verbose;
} LocusPocusOptions;

// Usage statement
void print_usage(FILE *outstream)
{
  fprintf( outstream,
"Usage: ./locuspocus [options] gff3file1 [gff3file2 gff3file3 ...]\n"
"  Options:\n"
"    -d|--debug            print detailed debugging messages to terminal\n"
"                          (standard error)\n"
"    -h|--help             print this help message and exit\n"
"    -n|--numprocs:INT     number of processors to utilize; default is 1\n"
"    -o|--outfile: FILE    name of file to which results will be written;\n"
"                          default is terminal (standard output)\n"
"    -v|--verbose          print detailed log messages to terminal (standard error)\n\n" );
}

// Main program
int main(int argc, char **argv)
{
  // Parse options from command line
  int opt = 0;
  int optindex = 0;
  const char *optstr = "dhn:o:v";
  const struct option locuspocus_options[] =
  {
    { "debug",    no_argument, NULL, 'd' },
    { "help",     no_argument, NULL, 'h' },
    { "numprocs", no_argument, NULL, 'n' },
    { "outfile",  no_argument, NULL, 'o' },
    { "verbose",  no_argument, NULL, 'v' },
  };
  LocusPocusOptions options = { 0, 0, stdout, 0 };
  for( opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex))
  {
    switch(opt)
    {
      case 'd':
        options.debug = 1;
        options.verbose = 1;
        break;
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'n':
        if(sscanf(optarg, "%d", &options.numprocs) == EOF)
        {
          fprintf(stderr, "[LocusPocus] error: could not convert number of "
                  "processors '%s' to an integer", optarg);
          exit(1);
        }
        break;
      case 'o':
        options.outstream = fopen(optarg, "w");
        if(options.outstream == NULL)
        {
          fprintf( stderr, "[LocusPocus] error: could not open outfile '%s'\n",
                   optarg );
          exit(1);
        }
        break;
      case 'v':
        options.verbose = 1;
        break;
    }
  }
  int numfiles = argc - optind;
  if(numfiles < 1)
  {
    fprintf( stderr, "[LocusPocus] error: must provide at least one GFF3 file "
             "as input\n");
    print_usage(stderr);
    exit(1);
  }

  // Load data into memory
  gt_lib_init();
  omp_set_num_threads(1);
  AgnLogger *logger = agn_logger_new();
  AgnLocusIndex *loci = agn_locus_index_new(false);
  unsigned long numloci = agn_locus_index_parse_disk(loci, numfiles,
                              (const char **)argv + optind, options.numprocs,
                              logger);
  if(options.verbose)
    fprintf(stderr, "[LocusPocus] found %lu total loci\n", numloci);
  bool haderror = agn_logger_print_all(logger, stderr, "[LocusPocus] loading "
                                       "features from %d input files",
                                       numfiles);
  if(haderror)
    return 1;
  agn_logger_delete(logger);

  // Iterate over each gene to find sets of mutually overlapping genes and print
  // them
  fputs("##gff-version\t3\n", options.outstream);
  GtStrArray *seqids = agn_locus_index_seqids(loci);
  if(options.verbose)
  {
    fprintf(stderr, "[LocusPocus] found %lu sequences\n",
            gt_str_array_size(seqids));
  }
  unsigned long i,j;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seqloci = agn_locus_index_get(loci, seqid);
    gt_array_sort(seqloci, (GtCompare)agn_locus_compare);
    if(options.verbose)
    {
      fprintf(stderr, "[LocusPocus] found %lu loci for sequence '%s'\n",
              gt_array_size(seqloci), seqid);
    }
    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnLocus *locus = *(AgnLocus **)gt_array_get(seqloci, j);
      agn_locus_print(locus, options.outstream, "AEGeAn::LocusPocus");
    }
    gt_array_delete(seqloci);
  }

  agn_locus_index_delete(loci);
  gt_lib_clean();
  return 0;
}
