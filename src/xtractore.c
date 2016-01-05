/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <string.h>
#include "genometools.h"
#include "aegean.h"

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

// Simple data structure for program options
typedef struct
{
  FILE *idfile;
  FILE *outfile;
  bool typeoverride;
  GtHashmap *typestoextract;
  bool verbose;
  unsigned width;
  bool debug;
} XtractoreOptions;

// Simple data structure to group genomic coordinates and strand together
typedef struct
{
  GtRange  r;
  GtStrand s;
} XtractRegion;


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Release memory held by program options/settings.
 */
static void xtract_options_free_memory(XtractoreOptions *options);

/**
 * @function Parse program settings from command line arugments.
 */
static void xtract_options_parse(int argc, char **argv,
                                 XtractoreOptions *options, GtError *error);

/**
 * @function Set default values for program options/settings.
 */
static void xtract_options_set_defaults(XtractoreOptions *options);

/**
 * @function Function for comparing genomic ranges, uses ``gt_region_compare``
 * under the hood.
 */
static int xtract_region_compare(XtractRegion *r1, XtractRegion *r2);

/**
 * @function Retrieve the subsequence of ``sequence`` corresponding to the
 * genomic feature encoded by ``gn``.
 */
static char *xt_extract_subsequence(GtGenomeNode *gn, const GtUchar *sequence,
                                    GtUword seqlength);

/**
 * @function Print sequence out to a file, ensuring each line of sequence is no
 * longer than the specified ``width``.
 */
static void xt_format_sequence(FILE *outstream, char *sequence, unsigned width);

/**
 * @function Retrieves the feature type, handling pseudonodes by returning the
 * type of its first child.
 */
static const char *xt_get_feature_type(GtFeatureNode *fn);

/**
 * @function Given an annotation of a possibly discontinuous genomic feature,
 * return an array containing each subrange individually.
 */
static GtArray *xt_get_regions(GtGenomeNode *gn);

/**
 * @function Given a feature encoded by ``gn``, print the sequence corresponding
 * to that feature.
 */
static void
xt_print_feature_sequence(GtGenomeNode *gn, const GtUchar *sequence,
                          GtUword seqlength, XtractoreOptions *options);

/**
 * @function Print the program's usage statement.
 */
static void xt_print_usage(FILE *outstream);


//------------------------------------------------------------------------------
// Function implementations
//------------------------------------------------------------------------------

static void xtract_options_free_memory(XtractoreOptions *options)
{
  if(options->idfile != NULL)
    fclose(options->idfile);
  fclose(options->outfile);
  gt_hashmap_delete(options->typestoextract);
}

static void xtract_options_parse(int argc, char **argv,
                                 XtractoreOptions *options, GtError *error)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "dhi:o:t:Vvw:";
  char *type;
  const struct option xtractore_options[] =
  {
    { "debug",    no_argument,       NULL, 'd' },
    { "help",     no_argument,       NULL, 'h' },
    { "idfile",   required_argument, NULL, 'i' },
    { "outfile",  required_argument, NULL, 'o' },
    { "type",     required_argument, NULL, 't' },
    { "verbose",  no_argument,       NULL, 'V' },
    { "version",  no_argument,       NULL, 'v' },
    { "width",    required_argument, NULL, 'w' },
    { NULL,       no_argument,       NULL,  0  },
  };
  for(opt = getopt_long(argc, argv + 0, optstr, xtractore_options, &optindex);
      opt != -1;
      opt = getopt_long(argc, argv + 0, optstr, xtractore_options, &optindex))
  {
    if(opt == 'd')
    {
      options->debug = true;
    }
    else if(opt == 'h')
    {
      xt_print_usage(stdout);
      exit(0);
    }
    else if(opt == 'i')
    {
      options->idfile = fopen(optarg, "r");
      if(options->idfile == NULL)
        gt_error_set(error, "could not open ID file '%s'", optarg);
    }
    else if(opt == 'o')
    {
      options->outfile = fopen(optarg, "w");
      if(options->outfile == NULL)
        gt_error_set(error, "could not open output file '%s'", optarg);
    }
    else if(opt == 't')
    {
      if(options->typeoverride == false)
      {
        gt_hashmap_delete(options->typestoextract);
        options->typestoextract = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                                 NULL);
        options->typeoverride = true;
      }
      type = gt_cstr_dup(optarg);
      gt_hashmap_add(options->typestoextract, type, type);
    }
    else if(opt == 'V')
    {
      options->verbose = true;
      break;
    }
    else if(opt == 'v')
    {
      agn_print_version("Xtractore", stdout);
      exit(0);
    }
    else if(opt == 'w')
    {
      if(sscanf(optarg, "%u", &options->width) == EOF)
      {
        gt_error_set(error, "could not convert width '%s' to an integer",
                     optarg);
      }
    }
  }
}

static void xtract_options_set_defaults(XtractoreOptions *options)
{
  options->idfile = NULL;
  options->outfile = stdout;
  options->typeoverride = false;
  options->typestoextract = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  char *defaulttype = gt_cstr_dup("gene");
  gt_hashmap_add(options->typestoextract, defaulttype, defaulttype);
  options->verbose = false;
  options->width = 80;
  options->debug = false;
}

static int xtract_region_compare(XtractRegion *r1, XtractRegion *r2)
{
  return gt_range_compare(&r1->r, &r2->r);
}

static char *xt_extract_subsequence(GtGenomeNode *gn, const GtUchar *sequence,
                                    GtUword seqlength)
{
  char *outseq, *outseqp;
  GtArray *regions;
  GtStrand fstrand = GT_STRAND_UNKNOWN;
  GtUword i, nregions, length;

  regions = xt_get_regions(gn);
  nregions = gt_array_size(regions);
  length = 0;
  for(i = 0; i < nregions; i++)
  {
    XtractRegion *region = gt_array_get(regions, i);
    length += gt_range_length(&region->r);
    if(i == 0)
    {
      fstrand = region->s;
      continue;
    }
    if(region->s != fstrand)
    {
      GtStr *seqid = gt_genome_node_get_seqid(gn);
      fprintf(stderr, "[xtractore] error: feature at %s[%lu, %lu] belongs to "
              "a multifeature with strand inconsistencies\n", gt_str_get(seqid),
              region->r.start, region->r.end);
      exit(1);
    }
    if(region->r.end > seqlength)
    {
      GtStr *seqid = gt_genome_node_get_seqid(gn);
      fprintf(stderr, "[xtractore] error: feature at %s[%lu, %lu] exceeds "
              "sequence length of %lu\n", gt_str_get(seqid), region->r.start,
              region->r.end, seqlength);
      exit(1);
    }
  }

  if(fstrand == GT_STRAND_REVERSE && gt_array_size(regions) > 1)
    gt_array_reverse(regions);

  outseq = gt_malloc( sizeof(char) * (length + 1) );
  // Added following line to suppress Valgrind warning about uninitialized
  // variable. It complains since the subsequent block operates only on
  // ``outseqp`` and not ``outseq``.
  outseq[length] = '\0';
  outseqp = outseq;
  GtError *error = gt_error_new();
  for(i = 0; i < nregions; i++)
  {
    XtractRegion *region = gt_array_get(regions, i);
    GtUword rlength = gt_range_length(&region->r);
    strncpy(outseqp, (char *)(sequence + region->r.start - 1), rlength);
    if(region->s == GT_STRAND_REVERSE)
      gt_reverse_complement(outseqp, rlength, error);
    outseqp += rlength;
  }
  gt_error_delete(error);

  gt_array_delete(regions);
  return outseq;
}

static void xt_format_sequence(FILE *outstream, char *sequence, unsigned width)
{
  GtUword i;
  GtUword seqlen = strlen(sequence);
  for(i = 0; i < seqlen; i++)
  {
    if(i > 0 && i % width == 0 && width != 0)
      fputc('\n', outstream);
    fputc(sequence[i], outstream);
  }
  fputc('\n', outstream);
}

static const char *xt_get_feature_type(GtFeatureNode *fn)
{
  if(gt_feature_node_is_pseudo(fn))
  {
    GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(fn);
    GtFeatureNode *child = gt_feature_node_iterator_next(iter);
    const char *type = gt_feature_node_get_type(child);
    gt_feature_node_iterator_delete(iter);
    return type;
  }
  return gt_feature_node_get_type(fn);
}

static GtArray *xt_get_regions(GtGenomeNode *gn)
{
  GtArray *regions;
  GtFeatureNode *fn;
  GtRange r;
  GtStrand s;
  XtractRegion region;

  regions = gt_array_new( sizeof(XtractRegion) );
  fn = gt_feature_node_cast(gn);

  if(gt_feature_node_is_pseudo(fn))
  {
    GtFeatureNodeIterator *iter;
    GtFeatureNode *current;

    iter = gt_feature_node_iterator_new_direct(fn);
    for(current  = gt_feature_node_iterator_next(iter);
        current != NULL;
        current  = gt_feature_node_iterator_next(iter))
    {
      r = gt_genome_node_get_range((GtGenomeNode *)current);
      s = gt_feature_node_get_strand(current);
      region.r = r;
      region.s = s;
      gt_array_add(regions, region);
    }
    gt_feature_node_iterator_delete(iter);
    gt_array_sort(regions, (GtCompare)xtract_region_compare);
    return regions;
  }

  r = gt_genome_node_get_range((GtGenomeNode *)fn);
  s = gt_feature_node_get_strand(fn);
  region.r = r;
  region.s = s;
  gt_array_add(regions, region);
  return regions;
}

static void
xt_print_feature_sequence(GtGenomeNode *gn, const GtUchar *sequence,
                          GtUword seqlength, XtractoreOptions *options)
{
  char subseqid[1024];

  GtFeatureNode *fn = gt_feature_node_cast(gn);
  GtRange range = gt_genome_node_get_range(gn);
  GtStr *seqid = gt_genome_node_get_seqid(gn);
  const char *type = xt_get_feature_type(fn);
  agn_assert(type);

  GtStrand strand = gt_feature_node_get_strand(fn);
  sprintf(subseqid, "%s_%lu-%lu%c", gt_str_get(seqid), range.start, range.end,
          GT_STRAND_CHARS[strand]);
  const char *featlabel = agn_feature_node_get_label(fn);
  if(gt_feature_node_is_pseudo(fn))
  {
    GtFeatureNodeIterator *it = gt_feature_node_iterator_new_direct(fn);
    GtFeatureNode *child = gt_feature_node_iterator_next(it);
    featlabel = agn_feature_node_get_label(child);
    gt_feature_node_iterator_delete(it);
  }
  fprintf(options->outfile, ">%s %s\n", featlabel, subseqid);

  char *feat_seq = xt_extract_subsequence(gn, sequence, seqlength);
  if(strcmp(type, "CDS") == 0 &&
     strncmp(feat_seq, "ATG", 3) != 0 &&
     options->verbose)
  {
    fprintf(stderr, "[xtractore] warning: CDS at '%s' does not begin with "
            "ATG\n", subseqid);
  }
  xt_format_sequence(options->outfile, feat_seq, options->width);
  gt_free(feat_seq);
}

static void xt_print_usage(FILE *outstream)
{
  fprintf(outstream,
"\nxtractore: extract sequences corresponding to annotated features from the\n"
"           given sequence file\n\n"
"Usage: xtractore [options] features.gff3 sequences.fasta\n"
"  Options:\n"
"    -d|--debug            print debugging output\n"
"    -h|--help             print this help message and exit\n"
"    -i|--idfile: FILE     file containing a list of feature IDs (1 per line\n"
"                          with no spaces); if provided, only features with\n"
"                          IDs in this file will be extracted\n"
"    -o|--outfile: FILE    file to which output sequences will be written;\n"
"                          default is terminal (stdout)\n"
"    -t|--type: STRING     feature type to extract; can be used multiple\n"
"                          times to extract features of multiple types\n"
"    -v|--version          print version number and exit\n"
"    -V|--verbose          print verbose warning and error messages\n"
"    -w|--width: INT       width of each line of sequence in the Fasta\n"
"                          output; default is 80; set to 0 for no\n"
"                          formatting\n\n");
}

int main(int argc, char **argv)
{
  const char *featfile, *seqfile;
  char *seqdesc;
  GtError *error;
  GtHashmap *seqs_observed;
  GtFeatureIndex *features;
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams;
  GtSeqIterator *seqiter;
  GtStrArray *seqfastas;
  const GtUchar *sequence;
  GtUword seqlength;
  int result;
  gt_lib_init();

  XtractoreOptions options;
  xtract_options_set_defaults(&options);
  error = gt_error_new();
  xtract_options_parse(argc, argv, &options, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "[xtractore] error: %s", gt_error_get(error));
    return 1;
  }
  int numfiles = argc - optind;
  if(numfiles < 2)
  {
    fprintf(stderr, "[xtractore] error: must provide a feature annotation file "
            "(GFF3 format) and a sequence file (Fasta format)\n");
    xt_print_usage(stderr);
    return 1;
  }
  featfile = argv[optind + 0];
  seqfile  = argv[optind + 1];

  streams = gt_queue_new();

  current_stream = gt_gff3_in_stream_new_unsorted(1, &featfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_filter_stream_new(last_stream, options.typestoextract);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  features = gt_feature_index_memory_new();
  current_stream = gt_feature_out_stream_new(last_stream, features);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[xtractore] error processing GFF3: %s\n",
            gt_error_get(error));
    return 1;
  }

  seqs_observed = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
  seqfastas = gt_str_array_new();
  gt_str_array_add_cstr(seqfastas, seqfile);
  seqiter = gt_seq_iterator_sequence_buffer_new(seqfastas, error);
  GtUword featcounter = 0;
  while((result = gt_seq_iterator_next(seqiter, &sequence, &seqlength, &seqdesc,
                                       error)) > 0)
  {
    char *seqid = strtok(seqdesc, " \n\t");
    gt_hashmap_add(seqs_observed, gt_cstr_dup(seqid), seqs_observed);
    GtArray *seqfeatures =
                gt_feature_index_get_features_for_seqid(features, seqid, error);
    GtUword nfeats = gt_array_size(seqfeatures);
    if(nfeats == 0)
    {
      gt_array_delete(seqfeatures);
      continue;
    }

    if(nfeats > 1)
      gt_array_sort(seqfeatures, (GtCompare)agn_genome_node_compare);

    GtUword i;
    for(i = 0; i < nfeats; i++)
    {
      GtGenomeNode *gn = *(GtGenomeNode **)gt_array_get(seqfeatures, i);
      xt_print_feature_sequence(gn, sequence, seqlength, &options);
      featcounter += 1;
      if(featcounter % 1000 == 0 && options.debug)
        fputs("..........", stderr);
      if(featcounter % 10000 == 0 && options.debug)
        fputs("\n", stderr);
    }
    gt_array_delete(seqfeatures);
  }
  if(result == -1)
  {
    fprintf(stderr, "[xtractore] error processing Fasta: %s\n",
            gt_error_get(error));
    return 1;
  }
  if(featcounter >= 1000 && options.debug)
    fputs("\n", stderr);

  GtStrArray *gff3seqids = gt_feature_index_get_seqids(features, error);
  GtUword i;
  for(i = 0; i < gt_str_array_size(gff3seqids); i++)
  {
    const char *seqid = gt_str_array_get(gff3seqids, i);
    if(gt_hashmap_get(seqs_observed, seqid) == NULL)
    {
      fprintf(stderr, "[AEGeAn::Xtractore] warning: sequence '%s' contains "
              "annotated features but no sequence was provided\n", seqid);
    }
  }
  gt_str_array_delete(gff3seqids);
  gt_hashmap_delete(seqs_observed);

  gt_seq_iterator_delete(seqiter);
  gt_str_array_delete(seqfastas);
  while(gt_queue_size(streams) > 0)
  {
    GtNodeStream *stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);
  gt_feature_index_delete(features);
  gt_error_delete(error);
  xtract_options_free_memory(&options);
  gt_lib_clean();
  return 0;
}
