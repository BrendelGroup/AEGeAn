#include <omp.h>
#include <string.h>
#include <time.h>
#include "AgnPairwiseCompareLocus.h"
#include "AgnUtils.h"
#include "PeNodeVisitor.h"
#include "PeReports.h"

// Main method
int main(int argc, char * const argv[])
{
  // Parse command-line arguments
  PeOptions options;
  pe_set_option_defaults(&options);
  int optind = pe_parse_options(argc, argv, &options);
  options.refrfile = argv[optind];
  options.predfile = argv[optind + 1];
  if(options.refrfile == NULL || options.predfile == NULL)
  {
    fprintf(stderr, "[ParsEval] error: could not parse input filenames\n");
    return EXIT_FAILURE;
  }

  // Grab start time
  time_t start_time;
  struct tm *start_time_info;
  char time_buffer[128];
  time(&start_time);
  start_time_info = localtime(&start_time);
  strftime(time_buffer, 128, "%d %b %Y, %I:%M%p", start_time_info);

  // Initialize ParsEval
  gt_lib_init();
  GtTimer *timer_all = gt_timer_new();
  gt_timer_start(timer_all);
  GtTimer *timer_short = gt_timer_new();
  fputs("[ParsEval] Begin ParsEval\n", stderr);


  //----- Load data into memory -----
  //---------------------------------
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin loading data\n", stderr);

  GtError *error = gt_error_new();
  AgnGeneValidator *validator = agn_gene_validator_new();
  GtGenomeNode *gn;
  bool loaderror;

  GtNodeStream *gff3in = gt_gff3_in_stream_new_unsorted(1, &options.refrfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtFeatureIndex *refrfeats = gt_feature_index_memory_new();
  GtNodeVisitor *refrvisitor = pe_node_visitor_new(refrfeats, validator);
  while(!(loaderror = gt_node_stream_next(gff3in, &gn, error)) && gn)
  {
    gt_genome_node_accept(gn, refrvisitor, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "%s\n", gt_error_get(error));
      return EXIT_FAILURE;
    }
  }
  gt_node_stream_delete(gff3in);
  gt_node_visitor_delete(refrvisitor);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "%s\n", gt_error_get(error));
    return EXIT_FAILURE;
  }

  gff3in = gt_gff3_in_stream_new_unsorted(1, &options.predfile);
  gt_gff3_in_stream_check_id_attributes((GtGFF3InStream *)gff3in);
  gt_gff3_in_stream_enable_tidy_mode((GtGFF3InStream *)gff3in);
  GtFeatureIndex *predfeats = gt_feature_index_memory_new();
  GtNodeVisitor *predvisitor = pe_node_visitor_new(predfeats, validator);
  while(!(loaderror = gt_node_stream_next(gff3in, &gn, error)) && gn)
  {
    gt_genome_node_accept(gn, predvisitor, error);
    if(gt_error_is_set(error))
    {
      fprintf(stderr, "%s\n", gt_error_get(error));
      return EXIT_FAILURE;
    }
  }
  gt_node_stream_delete(gff3in);
  gt_node_visitor_delete(predvisitor);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "%s\n", gt_error_get(error));
    return EXIT_FAILURE;
  }

  GtStrArray *refrseqids = gt_feature_index_get_seqids(refrfeats, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "%s\n", gt_error_get(error));
    return EXIT_FAILURE;
  }
  GtStrArray *predseqids = gt_feature_index_get_seqids(predfeats, error);
  if(gt_error_is_set(error))
  {
    fprintf(stderr, "%s\n", gt_error_get(error));
    return EXIT_FAILURE;
  }
  GtStrArray *seqids = agn_seq_intersection(refrseqids, predseqids);

  gt_str_array_delete(refrseqids);
  gt_str_array_delete(predseqids);
  gt_error_delete(error);
  gt_timer_stop(timer_short);
  gt_timer_show_formatted( timer_short, "[ParsEval] Finished loading data (%ld.%06ld seconds)\n",
                           stderr );
  unsigned long numseqs = gt_str_array_size(seqids);


  //----- Parse loci -----
  //----------------------
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin parsing loci\n", stderr);

  int i;
  unsigned long total = 0;
  int rank;
  GtArray *loci = gt_array_new( sizeof(GtArray *) );
  GtArray **p_list = (GtArray **)gt_malloc( numseqs * sizeof(GtArray *) );
  #pragma omp parallel private(i, rank)
  {
    rank = omp_get_thread_num();

    #pragma omp for schedule(static)
    for(i = 0; i < numseqs; i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      GtArray *seq_loci = agn_parse_loci(seqid, refrfeats, predfeats);
      #pragma omp critical
      {
        total += gt_array_size(seq_loci);
        p_list[i] = seq_loci;
        fprintf( stderr, "[ParsEval]     loci on seq '%s' identified by processor %d\n",
                 seqid, rank );
      }
    }
  } // End parallelize
  for(i = 0; i < numseqs; i++)
  {
    gt_array_add(loci, p_list[i]);
  }
  gt_free(p_list);
  gt_timer_stop(timer_short);
  gt_timer_show_formatted( timer_short, "[ParsEval] Finished parsing loci (%ld.%06ld seconds)\n",
                           stderr );
  gt_feature_index_delete(refrfeats);
  gt_feature_index_delete(predfeats);

  // Stop now if there are no loci
  if(total == 0)
  {
    fprintf(stderr, "[ParsEval] Warning: found no loci to analyze");
    gt_timer_stop(timer_all);
    gt_timer_show_formatted( timer_all,"[ParsEval] ParsEval complete! (total runtime: %ld.%06ld "
                             "seconds)\n\n", stderr );
    fclose(options.outfile);
    for(i = 0; i < gt_array_size(loci); i++)
    {
      GtArray *seq_loci = *(GtArray **)gt_array_get(loci, i);
      gt_array_delete(seq_loci);
    }
    gt_array_delete(loci);
    gt_str_array_delete(seqids);
    gt_timer_delete(timer_all);
    gt_timer_delete(timer_short);
    if(gt_lib_clean() != 0)
    {
      fputs("error: issue cleaning GenomeTools library\n", stderr);
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }


  //----- Locus comparison -----
  //----------------------------
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin comparative analysis\n", stderr);

  // Counts and stats needed to print summary report
  AgnSummaryData summary_data;
  agn_summary_data_init(&summary_data);

  // Counts and stats at the sequence level
  AgnSummaryData *seqlevel_summary_data = gt_malloc( sizeof(AgnSummaryData) * numseqs );
  for(i = 0; i < numseqs; i++)
  {
    agn_summary_data_init(&seqlevel_summary_data[i]);
  }

  if(strcmp(options.outfmt, "csv") == 0)
  {
    // Print header
    fputs( "Sequence,Start,End,"
           "Reference Transcript(s),Prediction Transcript(s),"
           "Reference CDS segments,Prediction CDS segments,"
           "Correct CDS segments,Missing CDS segments,Wrong CDS segments,"
           "CDS structure sensitivity,CDS structure specificity,"
           "CDS structure F1,CDS structure AED,"
           "Reference exons,Prediction exons,"
           "Correct exons,Missing exons,Wrong exons,"
           "Exon sensitivity,Exon specificity,"
           "Exon F1,Exon AED,"
           "Reference UTR segments,Prediction UTR segments,"
           "Correct UTR segments,Missing UTR segments,Wrong UTR segments,"
           "UTR structure sensitivity,UTR structure specificity,"
           "UTR structure F1,UTR structure AED,Overall identity,"
           "CDS nucleotide matching coefficient,CDS nucleotide correlation coefficient,"
           "CDS nucleotide sensitivity,CDS nucleotide specificity,"
           "CDS nucleotide F1,CDS nucleotide AED,"
           "UTR nucleotide matching coefficient,UTR nucleotide correlation coefficient,"
           "UTR nucleotide sensitivity,UTR nucleotide specificity,"
           "UTR nucleotide F1,UTR nucleotide AED\n", options.outfile );
  }

  // Iterate over all the sequences
  for(i = 0; i < numseqs; i++)
  {
    // For text output, write locus results to a temp file
    // For HTML output, each sequence and each locus gets a dedicated .html file
    FILE *seqfile = NULL;
    FILE *locusgff3file = NULL;
    GtArray *seq_loci = *(GtArray **)gt_array_get(loci, i);
    AgnPairwiseCompareLocusSummary *locus_summaries = (AgnPairwiseCompareLocusSummary *)
                                  gt_malloc( sizeof(AgnPairwiseCompareLocusSummary) * gt_array_size(seq_loci) );

    if(!options.summary_only)
    {
      const char *seqid = gt_str_array_get(seqids, i);

      char filename[512];
      if(options.html)
      {
        // Create output directory for this sequence
        char dircmd[512];
        sprintf(dircmd, "mkdir %s/%s", options.outfilename, seqid);
        if(system(dircmd) != 0)
        {
          fprintf(stderr, "error: could not open directory '%s/%s'\n", options.outfilename, seqid);
          exit(1);
        }
        if(options.debug)
          fprintf(stderr, "debug: opening directory '%s'\n", dircmd);

        // Create symbolic link to LICENSE in root directory
        sprintf( dircmd, "ln -s ../LICENSE %s/%s/LICENSE", options.outfilename,
                 seqid );
        if(system(dircmd) != 0)
        {
          fprintf(stderr, "warning: could not create symbolic link to LICENSE\n");
        }

        // Create summary page for this sequence
        sprintf(filename, "%s/%s/index.html", options.outfilename, seqid);
        if(options.debug)
          fprintf(stderr, "debug: opening outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w");
        pe_print_seqfile_header(seqfile, seqid);
      }
      else
      {
        sprintf(filename, "%s.%s", options.outfilename, seqid);
        if(options.debug)
          fprintf(stderr, "debug: opening temp outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w");
      }
    }

    if(options.locusgff3)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char filename[512];
      sprintf(filename, "%s.%s", options.locusfilename, seqid);
      if(options.debug)
        fprintf(stderr, "debug: opening temp locus file '%s'\n", filename);
      locusgff3file = agn_fopen(filename, "w");
    }

    // Distribute the loci for this sequence across the p processors
    int j;
    #pragma omp parallel private(rank, j)
    {
      AgnSummaryData summary_data_local;
      agn_summary_data_init(&summary_data_local);

      rank = omp_get_thread_num();

      #pragma omp for schedule(dynamic)
      for(j = 0; j < gt_array_size(seq_loci); j++)
      // Begin parallelize loop
      // Parallel loop over the loci, not the clique pairs
      {
        agn_comparison_counts_init(&locus_summaries[j].counts);
        AgnPairwiseCompareLocus *locus = *(AgnPairwiseCompareLocus **)gt_array_get(seq_loci, j);

        GtArray *clique_pairs = agn_pairwise_compare_locus_get_clique_pairs(locus, options.trans_per_locus);
        unsigned long comparisons = 0;
        if(clique_pairs != NULL)
          comparisons = gt_array_size(clique_pairs);
// fprintf(stderr, "DELETEME clique pairs: %lu\n", comparisons);
// fprintf(stderr, "DELETEME one\n");

        if(options.complimit != 0 && comparisons > options.complimit)
        {
          if(options.debug)
            fprintf( stderr, "debug: locus %s[%lu, %lu] contains %lu transcript (or transcript"
                     " clique) pairs, exceeds the limit of %d, moving on\n",
                     agn_pairwise_compare_locus_get_seqid(locus), agn_pairwise_compare_locus_get_start(locus),
                     agn_pairwise_compare_locus_get_end(locus), comparisons, options.complimit );
        }
        else if(agn_pairwise_compare_locus_filter(locus, &options.filters))
        {
          locus_summaries[j].start = 0;
          locus_summaries[j].end   = 0;
          if(options.debug)
            fprintf( stderr, "debug: locus %s[%lu, %lu] was filtered, moving on\n",
                     agn_pairwise_compare_locus_get_seqid(locus), agn_pairwise_compare_locus_get_start(locus),
                     agn_pairwise_compare_locus_get_end(locus));
        }
        else
        {
          summary_data_local.counts.num_loci++;

          // Grab some summary counts
          summary_data_local.counts.refr_genes       += agn_pairwise_compare_locus_num_refr_genes(locus);
          summary_data_local.counts.refr_transcripts += agn_pairwise_compare_locus_num_refr_transcripts(locus);
          summary_data_local.counts.pred_genes       += agn_pairwise_compare_locus_num_pred_genes(locus);
          summary_data_local.counts.pred_transcripts += agn_pairwise_compare_locus_num_pred_transcripts(locus);
          if(agn_pairwise_compare_locus_num_refr_genes(locus) == 0)
            summary_data_local.counts.unique_pred++;
          else if(agn_pairwise_compare_locus_num_pred_genes(locus) == 0)
            summary_data_local.counts.unique_refr++;

          if(options.html && !options.summary_only)
          {
            locus_summaries[j].start = agn_pairwise_compare_locus_get_start(locus);
            locus_summaries[j].end = agn_pairwise_compare_locus_get_end(locus);
            locus_summaries[j].length = agn_pairwise_compare_locus_get_length(locus);
            locus_summaries[j].refr_transcripts = agn_pairwise_compare_locus_num_refr_transcripts(locus);
            locus_summaries[j].pred_transcripts = agn_pairwise_compare_locus_num_pred_transcripts(locus);
            locus_summaries[j].total = comparisons;
          }

          if( clique_pairs != NULL &&
              (options.complimit == 0 || comparisons <= options.complimit) )
          {
            unsigned long k;
            GtTimer *comp_timer = gt_timer_new();
            gt_timer_start(comp_timer);
            for(k = 0; k < gt_array_size(clique_pairs); k++)
            {
// fprintf(stderr, "DELETEME k=%lu\n", k);
              AgnCliquePair *pair = *(AgnCliquePair **)gt_array_get(clique_pairs, k);
              agn_clique_pair_build_model_vectors(pair);
              agn_clique_pair_comparative_analysis(pair);
            }
            gt_timer_stop(comp_timer);

            // Begin maybe delete
            if(options.debug)
            {
              fprintf( stderr,
                       "debug: (proc %d) locus: %s[%lu, %lu]; length: %lu, trans: %lu,%lu; pairs: ",
                       rank, agn_pairwise_compare_locus_get_seqid(locus), agn_pairwise_compare_locus_get_start(locus),
                       agn_pairwise_compare_locus_get_end(locus), agn_pairwise_compare_locus_get_length(locus),
                       agn_pairwise_compare_locus_num_refr_transcripts(locus),
                       agn_pairwise_compare_locus_num_pred_transcripts(locus) );
              if(clique_pairs == NULL)
                fprintf(stderr, "0; no comparison");
              else if(options.complimit != 0 && gt_array_size(clique_pairs) > options.complimit)
                fprintf(stderr, "%lu; no comparison", gt_array_size(clique_pairs));
              else
              {
                fprintf(stderr, "%lu, ", gt_array_size(clique_pairs));
                gt_timer_show_formatted(comp_timer, "time: %ld.%06ld seconds", stderr);
              }
              fputs("\n", stderr);
            }
            // End maybe delete
            gt_timer_delete(comp_timer);

            if( clique_pairs != NULL &&
                (options.complimit == 0 || gt_array_size(clique_pairs) <= options.complimit) )
            {
              agn_pairwise_compare_locus_calc_splice_complexity_refr(locus);
              agn_pairwise_compare_locus_calc_splice_complexity_pred(locus);
            }

            AgnSummaryData data;
            agn_pairwise_compare_locus_get_summary_data(locus, &data);
            agn_pairwise_compare_locus_aggregate_results(locus, &data);
            agn_pairwise_compare_locus_aggregate_results(locus, &summary_data_local);
            GtArray *reported_pairs = agn_pairwise_compare_locus_find_best_pairs(locus);
            if(reported_pairs == NULL)
              locus_summaries[j].reported = 0;
            else
              locus_summaries[j].reported = gt_array_size(reported_pairs);
            locus_summaries[j].counts = data.counts;
          }

          if(!options.summary_only)
          {
            AgnSummaryData data;
            agn_pairwise_compare_locus_get_summary_data(locus, &data);
            if(options.html)
            {
              pe_gene_locus_print_results(locus, seqfile, &options);
            }
            else
            {
              #pragma omp critical(print_results)
              pe_gene_locus_print_results(locus, seqfile, &options);
            }
#ifndef WITHOUT_CAIRO
            if(options.locus_graphics)
            {
              AgnPairwiseCompareLocusPngMetadata metadata;
              pe_gene_locus_get_png_filename(locus, metadata.filename, options.outfilename);
              metadata.graphic_width = pe_gene_locus_get_graphic_width(locus);
              sprintf(metadata.stylefile, "%s/pe.style", options.data_path);
              metadata.refrfile = options.refrfile;
              metadata.predfile = options.predfile;
              metadata.refrlabel = options.refrlabel;
              metadata.predlabel = options.predlabel;
              metadata.track_order_func = pe_track_order;

              agn_pairwise_compare_locus_print_png(locus, &metadata);
            }
#endif
          }

          if(options.locusgff3)
            agn_pairwise_compare_locus_to_gff3(locus, locusgff3file);
        }
// fprintf(stderr, "DELETEME two\n");
        agn_pairwise_compare_locus_delete(locus);
// fprintf(stderr, "DELETEME three\n");
      } // End parallel for loop

      #pragma omp critical(aggregate_stats)
      {
        agn_summary_data_combine(&summary_data, &summary_data_local);
        agn_summary_data_combine(&seqlevel_summary_data[i], &summary_data_local);
      }
    } // End pragma omp parallel

    if(!options.summary_only)
    {
      if(options.html)
      {
        for(j = 0; j < gt_array_size(seq_loci); j++)
        {
          AgnPairwiseCompareLocusSummary *locus_summary = locus_summaries + j;
          if(locus_summary->start > 0 && locus_summary->end > 0)
            pe_print_locus_to_seqfile( seqfile, locus_summary->start, locus_summary->end,
                                       locus_summary->length, locus_summary->refr_transcripts,
                                       locus_summary->pred_transcripts, &locus_summary->counts );
        }
        pe_print_seqfile_footer(seqfile);
      }
      fclose(seqfile);
    }

    if(options.locusgff3)
      fclose(locusgff3file);
    gt_free(locus_summaries);
    gt_array_delete(seq_loci);
  } // End iterate over sequences

  gt_timer_stop(timer_short);
  gt_timer_show_formatted( timer_short, "[ParsEval] Finished comparative analysis (%ld.%06ld "
                           "seconds)\n", stderr );


  //----- Print summary statistics, combine output -----
  //----------------------------------------------------
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin printing summary, combining output\n", stderr);

  pe_print_summary( time_buffer, argc, argv, seqids, &summary_data, seqlevel_summary_data,
                    options.outfile, &options );
  if(options.outfile == stdout)
    fflush(stdout);
  else
    fclose(options.outfile);
  if(!options.html && !options.summary_only)
  {
    int result;
    for(i = 0; i < numseqs; i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char command[1024];
      sprintf( command, "cat %s.%s >> %s && rm %s.%s", options.outfilename, seqid,
               options.outfilename, options.outfilename, seqid );
      if(options.outfile == stdout)
      {
        sprintf( command, "cat %s.%s && rm %s.%s", options.outfilename, seqid, options.outfilename,
                 seqid );
      }
      if(options.debug)
        fprintf(stderr, "debug: merging output files: %s\n", command);
      result = system(command);
      if(result)
      {
        fprintf(stderr, "[ParsEval] error: issue merging GFF3 files: %s\n", command);
        exit(1);
      }
    }
  }
  if(options.locusgff3)
  {
    int result;
    FILE *gff3out = agn_fopen(options.locusfilename, "w");
    fputs
    (
      "##gff-version\t3\n"
      "#\n"
      "# If visualizing with GBrowse, add the following to your track configuration\n"
      "# and customize as desired. See the ParsEval README for more details.\n"
      "#\n"
      "# glyph       = heat_map\n"
      "# min_score   = 0\n"
      "# max_score   = 1\n"
      "# start_color = red\n"
      "# end_color   = green\n"
      "# link        = http://localhost/myparsevalhtml/$ref/$start-$end.html\n"
      "#\n",
      gff3out
    );
    fclose(gff3out);

    for(i = 0; i < numseqs; i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char command[1024];
      sprintf( command, "cat %s.%s >> %s && rm %s.%s", options.locusfilename, seqid,
               options.locusfilename, options.locusfilename, seqid );
      if(options.debug)
        fprintf(stderr, "debug: merging GFF3 files: %s\n", command);
      result = system(command);
      if(result)
      {
        fprintf(stderr, "[ParsEval] error: issue merging GFF3 files: %s\n", command);
        exit(1);
      }
    }
  }
  gt_timer_stop(timer_short);
  gt_timer_show_formatted( timer_short, "[ParsEval] Finished printing summary, combining output "
                           "(%ld.%06ld seconds)\n", stderr );
  if(options.outfile == stdout)
    fclose(stdout);


  // All done!
  gt_timer_stop(timer_all);
  gt_timer_show_formatted( timer_all,"[ParsEval] ParsEval complete! (total runtime: %ld.%06ld "
                           "seconds)\n\n", stderr );

  // Free up remaining memory
  gt_free(seqlevel_summary_data);
  gt_array_delete(loci);
  gt_str_array_delete(seqids);
  gt_timer_delete(timer_all);
  gt_timer_delete(timer_short);
  if(gt_lib_clean() != 0)
  {
    fputs("error: issue cleaning GenomeTools library\n", stderr);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
