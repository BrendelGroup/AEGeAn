#include <getopt.h>
#include "genometools.h"
#include "AgnGeneLocus.h"
#include "AgnLocusIndex.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  bool intloci;
  unsigned long delta;
  int numprocs;
  FILE *outstream;
  bool skipends;
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
"    -i|--intloci          parse and report interval loci rather than gene\n"
"                          loci\n"
"    -l|--delta: INT       when parsing interval loci, use the following\n"
"                          delta to extend gene loci and include potential\n"
"                          regulatory regions; default is 500\n"
"    -n|--numprocs: INT    number of processors to utilize; default is 1\n"
"    -o|--outfile: FILE    name of file to which results will be written;\n"
"                          default is terminal (standard output)\n"
"    -s|--skipends         when enumerating interval loci, exclude gene-less\n"
"                          iloci at either end of the sequence\n"
"    -v|--verbose          print detailed log messages to terminal (standard\n"
"                          error)\n\n" );
}

// Main program
int main(int argc, char **argv)
{
  // Parse options from command line
  int opt = 0;
  int optindex = 0;
  const char *optstr = "dhil:n:o:sv";
  const struct option locuspocus_options[] =
  {
    { "debug",    no_argument,       NULL, 'd' },
    { "help",     no_argument,       NULL, 'h' },
    { "intloci",  no_argument,       NULL, 'i' },
    { "delta",    required_argument, NULL, 'l' },
    { "numprocs", required_argument, NULL, 'n' },
    { "outfile",  required_argument, NULL, 'o' },
    { "skipends", no_argument,       NULL, 's' },
    { "verbose",  no_argument,       NULL, 'v' },
  };
  LocusPocusOptions options = { 0, 0, 500, 0, stdout, 0, 0 };
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
      case 'i':
        options.intloci = 1;
        break;
      case 'l':
        if(sscanf(optarg, "%lu", &options.delta) == EOF)
        {
          fprintf(stderr, "[LocusPocus] error: could not convert delta '%s' to "
                  "an integer", optarg);
          exit(1);
        }
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
      case 's':
        options.skipends = 1;
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
  AgnLogger *logger = agn_logger_new();
  AgnLocusIndex *loci = agn_locus_index_new(true);
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
  unsigned long i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seqloci;
    if(options.intloci)
    {
       seqloci = agn_locus_index_interval_loci(loci, seqid, options.delta,
                                               options.skipends);
    }
    else
    {
      seqloci = agn_locus_index_get(loci, seqid);
    }
    gt_array_sort(seqloci, (GtCompare)agn_gene_locus_array_compare);
    gt_array_reverse(seqloci);
    if(options.verbose)
    {
      fprintf(stderr,"[LocusPocus] found %lu loci for sequence '%s'\n",
              gt_array_size(seqloci), seqid);
    }
    while(gt_array_size(seqloci) > 0)
    {
      AgnGeneLocus **locus = gt_array_pop(seqloci);
      agn_gene_locus_to_gff3(*locus, options.outstream, "AEGeAn::LocusPocus");
      if(options.intloci)
        agn_gene_locus_delete(*locus);
    }
    gt_array_delete(seqloci);
  }

  agn_locus_index_delete(loci);
  gt_lib_clean();
  return 0;
}
