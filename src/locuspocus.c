#include <getopt.h>
#include "genometools.h"
#include "AgnGeneLocus.h"
#include "AgnLocusIndex.h"

// Simple data structure for program options
typedef struct
{
  bool debug;
  FILE *genestream;
  bool intloci;
  unsigned long delta;
  FILE *outstream;
  bool skipends;
  FILE *transstream;
  bool verbose;
} LocusPocusOptions;

// Usage statement
void print_usage(FILE *outstream)
{
  fprintf( outstream,
"Usage: ./locuspocus [options] gff3file1 [gff3file2 gff3file3 ...]\n"
"  Options:\n"
"    -d|--debug             print detailed debugging messages to terminal\n"
"                           (standard error)\n"
"    -g|--genemap: FILE     print a mapping from each gene annotation to its\n"
"                           corresponding locus to the given file\n"
"    -h|--help              print this help message and exit\n"
"    -i|--intloci           parse and report interval loci rather than gene\n"
"                           loci\n"
"    -l|--delta: INT        when parsing interval loci, use the following\n"
"                           delta to extend gene loci and include potential\n"
"                           regulatory regions; default is 500\n"
"    -o|--outfile: FILE     name of file to which results will be written;\n"
"                           default is terminal (standard output)\n"
"    -s|--skipends          when enumerating interval loci, exclude gene-less\n"
"                           iloci at either end of the sequence\n"
"    -t|--transmap: FILE    print a mapping from each transcript annotation\n"
"                           to its corresponding locus to the given file\n"
"    -v|--verbose           print detailed log messages to terminal (standard\n"
"                           error)\n\n" );
}

void parse_options(int argc, char **argv, LocusPocusOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "dg:hil:n:o:st:v";
  const struct option locuspocus_options[] =
  {
    { "debug",    no_argument,       NULL, 'd' },
    { "genemap",  required_argument, NULL, 'g' },
    { "help",     no_argument,       NULL, 'h' },
    { "intloci",  no_argument,       NULL, 'i' },
    { "delta",    required_argument, NULL, 'l' },
    { "outfile",  required_argument, NULL, 'o' },
    { "skipends", no_argument,       NULL, 's' },
    { "transmap", required_argument, NULL, 't' },
    { "verbose",  no_argument,       NULL, 'v' },
  };
  for( opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex);
       opt != -1;
       opt = getopt_long(argc, argv + 0, optstr, locuspocus_options, &optindex))
  {
    switch(opt)
    {
      case 'd':
        options->debug = 1;
        options->verbose = 1;
        break;
      case 'g':
        options->genestream = fopen(optarg, "w");
        if(options->genestream == NULL)
        {
          fprintf(stderr, "[LocusPocus] error: could not open genemap file "
                  "'%s'\n", optarg);
          exit(1);
        }
        break;
      case 'h':
        print_usage(stdout);
        exit(0);
        break;
      case 'i':
        options->intloci = 1;
        break;
      case 'l':
        if(sscanf(optarg, "%lu", &options->delta) == EOF)
        {
          fprintf(stderr, "[LocusPocus] error: could not convert delta '%s' to "
                  "an integer", optarg);
          exit(1);
        }
        break;
      case 'o':
        options->outstream = fopen(optarg, "w");
        if(options->outstream == NULL)
        {
          fprintf( stderr, "[LocusPocus] error: could not open outfile '%s'\n",
                   optarg );
          exit(1);
        }
        break;
      case 's':
        options->skipends = 1;
        break;
      case 't':
        options->transstream = fopen(optarg, "w");
        if(options->transstream == NULL)
        {
          fprintf(stderr, "[LocusPocus] error: could not open transmap file "
                  "'%s'\n", optarg);
          exit(1);
        }
        break;
      case 'v':
        options->verbose = 1;
        break;
    }
  }
}

// Main program
int main(int argc, char **argv)
{
  // Parse options from command line
  LocusPocusOptions options = { 0, NULL, 0, 500, stdout, 0, NULL, 0 };
  parse_options(argc, argv, &options);
  int numfiles = argc - optind;
  if(numfiles < 1)
  {
    fprintf( stderr, "[LocusPocus] error: must provide at least one GFF3 file "
             "as input\n");
    print_usage(stderr);
    return 1;
  }

  // Load data into memory
  gt_lib_init();
  AgnLogger *logger = agn_logger_new();
  AgnLocusIndex *loci = agn_locus_index_new(true);
  unsigned long numloci = agn_locus_index_parse_disk(loci, numfiles,
                              (const char **)argv + optind, logger);
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
      if(options.genestream != NULL)
        agn_gene_locus_print_gene_mapping(*locus, options.genestream);
      if(options.transstream != NULL)
        agn_gene_locus_print_transcript_mapping(*locus, options.transstream);
      if(options.intloci)
        agn_gene_locus_delete(*locus);
    }
    gt_array_delete(seqloci);
  }

  fclose(options.outstream);
  if(options.genestream != NULL)  fclose(options.genestream);
  if(options.transstream != NULL) fclose(options.transstream);
  agn_locus_index_delete(loci);
  gt_lib_clean();
  return 0;
}
