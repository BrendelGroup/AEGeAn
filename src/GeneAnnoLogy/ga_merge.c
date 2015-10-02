/**

Copyright (c) 2010-2015, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include <getopt.h>
#include <string.h>
#include "ga_commands.h"
#include "genometools.h"
#include "aegean.h"

typedef enum
{
  GA_MERGE_COMP,
  GA_MERGE_REPL,
  GA_MERGE_DEFAULT
} GaMergeStrategy;

typedef struct
{
  const char *source;
  const char *repo;
  const char **filenames;
  GaMergeStrategy strategy;
  int numfiles;
} GaMergeOptions;

static void ga_merge_print_usage(FILE *outstream)
{
  fputs("\n"
"[GeneAnnoLogy::merge] merge a new annotation into the repo\n\n"
"Usage: geneannology merge [options] repo new-annot.gff3\n"
"  Options:\n"
"    -h|--help                      print this help message and exit\n"
"    -m|--merge-strategy: STRING    'complement' indicates that gene models\n"
"                                   from 'new-annot.gff3' that overlap with\n"
"                                   gene models in the repo should be\n"
"                                   ignored; 'replace' indicates that gene\n"
"                                   models from 'new-annot.gff3' that overlap\n"
"                                   with gene models in the repo should\n"
"                                   replace those in the repo\n"
"    -s|--source: STRING            replace the 'source' label (GFF3 2nd\n"
"                                   column) of the input with the given value\n"
"    -v|--version                   print version number and exit\n\n",
  outstream);
}

static void ga_merge_parse_options(int argc, char * const *argv,
                                  GaMergeOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hm:s:v";
  const struct option merge_options[] =
  {
    { "help",           no_argument,       NULL, 'h' },
    { "merge-strategy", required_argument, NULL, 'm' },
    { "source",         required_argument, NULL, 's' },
    { "version",        no_argument,       NULL, 'v' },
  };

  options->source = NULL;
  options->repo = NULL;
  options->strategy = GA_MERGE_DEFAULT;
  options->filenames = NULL;
  options->numfiles = 0;

  for(opt  = getopt_long(argc, argv, optstr, merge_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv, optstr, merge_options, &optindex))
  {
    if(opt == 'h')
    {
      ga_merge_print_usage(stdout);
      exit(0);
    }
    else if(opt == 'm')
    {
      if(strcmp(optarg, "complement") == 0)
        options->strategy = GA_MERGE_COMP;
      else if(strcmp(optarg, "replace") == 0)
        options->strategy = GA_MERGE_REPL;
      else
      {
        ga_merge_print_usage(stderr);
        fprintf(stderr, "error: unknown merge strategy '%s'\n", optarg);
        exit(1);
      }
    }
    else if(opt == 's')
    {
      options->source = optarg;
    }
    else if(opt == 'v')
    {
      agn_print_version("GeneAnnoLogy::merge", stdout);
      exit(0);
    }
    else
    {
      ga_merge_print_usage(stderr);
      fprintf(stderr, "error: unknown option '%c'", opt);
      exit(1);
    }
  }

  if(options->strategy == GA_MERGE_DEFAULT)
  {
    ga_merge_print_usage(stderr);
    fprintf(stderr, "error: please specify merge strategy\n");
    exit(1);
  }

  int numargs = argc - optind;
  if(numargs == 0)
    return;

  options->repo = argv[optind];
  options->numfiles = numargs - 1;
  if(numargs == 1)
    return;

  options->filenames = (const char **)(argv + optind + 1);
}

static void ga_merge_annotations(GtFeatureIndex *dest, GtFeatureIndex *src,
                                 GtError *error)
{
  agn_assert(dest && src && error);
  GtStrArray *seqids = gt_feature_index_get_seqids(src, error);
  GtUword i,j;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *feats = gt_feature_index_get_features_for_seqid(src, seqid, error);
    for(j = 0; j < gt_array_size(feats); j++)
    {
      GtGenomeNode **feature = gt_array_get(feats, j);
      GtRange frange = gt_genome_node_get_range(*feature);
      GtArray *overlap = gt_array_new( sizeof(GtFeatureNode *) );
      gt_feature_index_get_features_for_range(dest,overlap,seqid,&frange,error);
      if(gt_array_size(overlap) == 0)
      {
        GtFeatureNode *fn = *(GtFeatureNode **)feature;
        gt_feature_index_add_feature_node(dest, fn, error);
      }
      gt_array_delete(overlap);
    }
    gt_array_delete(feats);
  }
  gt_str_array_delete(seqids);
}

int ga_merge(int argc, char * const *argv)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  agn_assert(argc > 0);
  GaMergeOptions options;
  ga_merge_parse_options(argc, argv, &options);
  if(options.repo == NULL)
  {
    ga_merge_print_usage(stderr);
    return 0;
  }
  else if(options.numfiles == 0)
  {
    fprintf(stderr, "[GeneAnnoLogy::merge] warning: repo=%s, reading input from"
            " stdin\n", options.repo);
    current_stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  }
  else
  {
    current_stream = gt_gff3_in_stream_new_unsorted(options.numfiles,
                                                    options.filenames);
  }
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  if(options.source)
  {
    GtStr *src = gt_str_new_cstr(options.source);
    GtNodeVisitor *nv = gt_set_source_visitor_new(src);
    current_stream = gt_visitor_stream_new(last_stream, nv);
    gt_queue_add(streams, current_stream);
    last_stream = current_stream;
    gt_str_delete(src);
  }

  GtFeatureIndex *new_features = gt_feature_index_memory_new();
  current_stream = gt_feature_out_stream_new(last_stream, new_features);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy::merge] error processing node stream: %s\n",
            gt_error_get(error));
  }
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }



  GtStrArray *input_seqids = gt_feature_index_get_seqids(new_features, error);
  GtStrArray *repo_files = ga_repo_get_filenames_for_seqids(options.repo,
                                                            input_seqids,
                                                            error);
  GtUword numfiles = gt_str_array_size(repo_files);
  char **filenames = gt_malloc( numfiles * sizeof(char *) );
  GtUword i;
  for(i = 0; i < numfiles; i++)
    filenames[i] = (char *)gt_str_array_get(repo_files, i);
  current_stream =
       gt_gff3_in_stream_new_unsorted(numfiles, (const char **)filenames);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  GtFeatureIndex *features = gt_feature_index_memory_new();
  current_stream = gt_feature_out_stream_new(last_stream, features);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy::merge] error processing node stream: %s\n",
            gt_error_get(error));
  }
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_str_array_delete(repo_files);
  gt_free(filenames);



  agn_assert(options.strategy == GA_MERGE_COMP ||
             options.strategy == GA_MERGE_REPL);
  GtFeatureIndex *mergedfeats;
  if(options.strategy == GA_MERGE_COMP)
  {
    ga_merge_annotations(features, new_features, error);
    mergedfeats = features;
  }
  else
  {
    ga_merge_annotations(new_features, features, error);
    mergedfeats = new_features;
  }



  current_stream = gt_feature_in_stream_new(mergedfeats);
  gt_feature_in_stream_use_orig_ranges((GtFeatureInStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  agn_locus_stream_set_source((AgnLocusStream *)current_stream, "GeneAnnoLogy");
  agn_locus_stream_skip_iiLoci((AgnLocusStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_repo_stream_open(last_stream, options.repo, error);
  if(current_stream == NULL)
  {
    fprintf(stderr, "[GeneAnnoLogy::merge] error setting up repo: %s\n",
            gt_error_get(error));
    return -1;
  }
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy::merge] error processing node stream: %s\n",
            gt_error_get(error));
  }
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }

  gt_feature_index_delete(features);
  gt_feature_index_delete(new_features);
  gt_error_delete(error);
  gt_logger_delete(logger);
  while(gt_queue_size(streams) > 0)
  {
    current_stream = gt_queue_get(streams);
    gt_node_stream_delete(current_stream);
  }
  gt_queue_delete(streams);
  return 0;
}
