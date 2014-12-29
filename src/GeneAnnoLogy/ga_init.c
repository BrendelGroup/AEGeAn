/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/

#include "ga_commands.h"
#include "genometools.h"
#include "aegean.h"

void ga_init_print_usage(FILE *outstream)
{
  fputs("Usage: geneannology init [options] repo gff3file\n", outstream);
}

int ga_init(int argc, char * const *argv)
{
  GtNodeStream *current_stream, *last_stream;
  GtQueue *streams = gt_queue_new();
  GtLogger *logger = gt_logger_new(true, "", stderr);
  GtError *error = gt_error_new();

  if(argc == 0)
  {
    ga_init_print_usage(stderr);
    return -1;
  }
  else if(argc == 1)
  {
    current_stream = gt_gff3_in_stream_new_unsorted(0, NULL);
  }
  else
  {
    const char **filenames = (const char **)argv + 1;
    current_stream = gt_gff3_in_stream_new_unsorted(argc - 1, filenames);
  }
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)current_stream);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)current_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_gene_stream_new(last_stream, logger);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_locus_stream_new(last_stream, 0);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_node_delete_stream_new(last_stream);
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;

  current_stream = agn_repo_stream_new(last_stream, argv[0], error);
  if(current_stream == NULL)
  {
    fprintf(stderr, "[GeneAnnoLogy] error setting up repo: %s\n",
            gt_error_get(error));
    return -1;
  }
  gt_queue_add(streams, current_stream);
  last_stream = current_stream;
  
  int result = gt_node_stream_pull(last_stream, error);
  if(result == -1)
  {
    fprintf(stderr, "[GeneAnnoLogy] error processing node stream: %s\n",
            gt_error_get(error));
  }
  
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
