/**

Copyright (c) 2010-2016, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#include <getopt.h>
#include "genometools.h"
#include "AgnMrnaRepVisitor.h"
#include "AgnPseudogeneFixVisitor.h"

typedef struct
{
  bool infer_introns;
  FILE *outstream;
  bool fix_pseudogenes;
  bool locus_parent;
  FILE *mapstream;
} PmrnaOptions;

static void print_usage(FILE *outstream)
{
  fprintf(outstream,
"\npmrna: filter out all but the primary isoform from each gene of the input\n"
"Usage: pmrna [options] < annot.gff3 > new.gff3\n"
"  Options:\n"
"    -h|--help           print this help message and exit\n"
"    -i|--introns        flag indicating that introns are declared explicitly\n"
"                        and do not need to be inferred from exon features;\n"
"                        default is to infer introns\n"
"    -l|--locus          report a single representative mRNA for each locus\n"
"                        instead of each gene\n"
"    -m|--map: FILE      write each gene/mRNA mapping to the specified file\n"
"    -p|--pseudogenes    disable pseudogene detection and correction\n\n");
}

static void parse_options(int argc, char **argv, PmrnaOptions *options)
{
  int opt = 0;
  int optindex = 0;
  const char *optstr = "hilm:p";
  const struct option pmrna_options[] =
  {
    { "help",        no_argument,       NULL, 'h' },
    { "introns",     no_argument,       NULL, 'i' },
    { "locus",       no_argument,       NULL, 'l' },
    { "map",         required_argument, NULL, 'm' },
    { "pseudogenes", no_argument,       NULL, 'o' },
    { NULL,          no_argument,       NULL,  0  },
  };
  for(opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex);
      opt != -1;
      opt  = getopt_long(argc, argv + 0, optstr, pmrna_options, &optindex))
  {
    if(opt == 'h')
    {
      print_usage(stdout);
      exit(0);
    }
    else if(opt == 'i')
      options->infer_introns = false;
    else if(opt == 'l')
      options->locus_parent = true;
    else if(opt == 'm')
    {
      options->mapstream = fopen(optarg, "w");
      if(options->mapstream == NULL)
      {
        fprintf(stderr, "error: could not create mapfile '%s'", optarg);
        exit(1);
      }
    }
    else if(opt == 'p')
      options->fix_pseudogenes = false;
  }
}

int main(int argc, char **argv)
{
  GtError *error;
  GtNodeStream *stream, *last_stream;
  GtQueue *streams;
  PmrnaOptions options = { true, NULL, true, false, NULL };
  parse_options(argc, argv, &options);

  //----------
  // Set up the processing stream
  //----------
  gt_lib_init();
  streams = gt_queue_new();

  stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  if(options.fix_pseudogenes)
  {
    stream = agn_pseudogene_fix_stream_new(last_stream);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  if(options.infer_introns)
  {
    stream = gt_add_introns_stream_new(last_stream);
    gt_queue_add(streams, stream);
    last_stream = stream;
  }

  GtNodeVisitor *nv = agn_mrna_rep_visitor_new(options.mapstream);
  if(options.locus_parent)
  {
    agn_mrna_rep_visitor_set_parent_type((AgnMrnaRepVisitor *)nv, "locus");
  }
  if(options.mapstream)
  {
    if(options.locus_parent)
    {
      fprintf(options.mapstream, "piLocusID\tMrnaID\n");
    }
    else
    {
      fprintf(options.mapstream, "GeneID\tMrnaID\n");
    }
  }
  stream = gt_visitor_stream_new(last_stream, nv);
  gt_queue_add(streams, stream);
  last_stream = stream;

  stream = gt_gff3_out_stream_new(last_stream, NULL);
  gt_gff3_out_stream_retain_id_attributes((GtGFF3OutStream *)stream);
  gt_queue_add(streams, stream);
  last_stream = stream;

  //----------
  // Execute the processing stream
  //----------
  error = gt_error_new();
  int had_err = gt_node_stream_pull(last_stream, error);
  if(had_err)
    fprintf(stderr, "Error processing node stream: %s\n", gt_error_get(error));
  gt_error_delete(error);

  //----------
  // Free memory
  //----------
  while(gt_queue_size(streams) > 0)
  {
    stream = gt_queue_get(streams);
    gt_node_stream_delete(stream);
  }
  gt_queue_delete(streams);
  if(options.mapstream != NULL)
    fclose(options.mapstream);
  gt_lib_clean();
  return had_err;
}
