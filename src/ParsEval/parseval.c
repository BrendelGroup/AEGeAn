#include <omp.h>
#include <string.h>
#include <time.h>
#include "AgnCanonNodeVisitor.h"
#include "AgnLocusIndex.h"
#include "AgnGeneLocus.h"
#include "AgnUtils.h"
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

  // Load data into memory and parse loci
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin loading data and parsing loci\n", stderr);

  AgnLogger *logger = agn_logger_new();
  AgnLocusIndex *locusindex = agn_locus_index_new();
  unsigned long total = agn_locus_index_parse_pairwise_disk(locusindex,
                            options.refrfile, options.predfile,
                            options.numprocs, logger);
  bool haderror = agn_logger_print_all(logger, stderr, "[ParsEval] parsing "
                                       "annotations from files '%s' and '%s'",
                                       options.refrfile, options.predfile);
  if(haderror) return EXIT_FAILURE;
  agn_logger_unset(logger);
  
  GtStrArray *seqids = agn_locus_index_seqids(locusindex);
  unsigned long numseqs = gt_str_array_size(seqids);
  GtArray *loci = gt_array_new( sizeof(GtArray *) );
  int i;
  for(i = 0; i < numseqs; i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seq_loci = agn_locus_index_get(locusindex, seqid);
    gt_array_sort(seq_loci,(GtCompare)agn_gene_locus_array_compare);
    gt_array_add(loci, seq_loci);
  }
  gt_timer_stop(timer_short);
  gt_timer_show_formatted(timer_short, "[ParsEval] Finished loading data and "
                          "parsing loci (%ld.%06ld seconds)\n", stderr);

  // Stop now if there are no loci
  if(total == 0)
  {
    fprintf(stderr, "[ParsEval] Warning: found no loci to analyze");
    gt_timer_stop(timer_all);
    gt_timer_show_formatted(timer_all,"[ParsEval] ParsEval complete! (total "
                            "runtime: %ld.%06ld seconds)\n\n", stderr );
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
  omp_set_num_threads(options.numprocs); // FIXME
  gt_timer_start(timer_short);
  fputs("[ParsEval] Begin comparative analysis\n", stderr);

  // Counts and stats needed to print summary report
  PeCompEvaluation summary_data;
  pe_comp_evalutation_init(&summary_data);

  // Counts and stats at the sequence level
  PeCompEvaluation *seqlevel_summary_data = gt_malloc( sizeof(PeCompEvaluation) * numseqs );
  for(i = 0; i < numseqs; i++)
  {
    pe_comp_evalutation_init(&seqlevel_summary_data[i]);
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
    AgnGeneLocusSummary *locus_summaries = (AgnGeneLocusSummary *)
                                  gt_malloc( sizeof(AgnGeneLocusSummary) * gt_array_size(seq_loci) );

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
    int rank, j;
    #pragma omp parallel private(rank, j)
    {
      PeCompEvaluation summary_data_local;
      pe_comp_evalutation_init(&summary_data_local);

      rank = omp_get_thread_num();

      #pragma omp for schedule(dynamic)
      for(j = 0; j < gt_array_size(seq_loci); j++)
      // Begin parallelize loop
      // Parallel loop over the loci, not the clique pairs
      {
        agn_comp_summary_init(&locus_summaries[j].counts);
        AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seq_loci, j);

        GtArray *clique_pairs = agn_gene_locus_get_clique_pairs(locus, options.trans_per_locus);
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
                     agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
                     agn_gene_locus_get_end(locus), comparisons, options.complimit );
        }
        else if(agn_gene_locus_filter(locus, &options.filters))
        {
          locus_summaries[j].start = 0;
          locus_summaries[j].end   = 0;
          if(options.debug)
            fprintf( stderr, "debug: locus %s[%lu, %lu] was filtered, moving on\n",
                     agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
                     agn_gene_locus_get_end(locus));
        }
        else
        {
          summary_data_local.counts.num_loci++;

          // Grab some summary counts
          summary_data_local.counts.refr_genes       += agn_gene_locus_num_refr_genes(locus);
          summary_data_local.counts.refr_transcripts += agn_gene_locus_num_refr_transcripts(locus);
          summary_data_local.counts.pred_genes       += agn_gene_locus_num_pred_genes(locus);
          summary_data_local.counts.pred_transcripts += agn_gene_locus_num_pred_transcripts(locus);
          if(agn_gene_locus_num_refr_genes(locus) == 0)
            summary_data_local.counts.unique_pred++;
          else if(agn_gene_locus_num_pred_genes(locus) == 0)
            summary_data_local.counts.unique_refr++;

          if(options.html && !options.summary_only)
          {
            locus_summaries[j].start = agn_gene_locus_get_start(locus);
            locus_summaries[j].end = agn_gene_locus_get_end(locus);
            locus_summaries[j].length = agn_gene_locus_get_length(locus);
            locus_summaries[j].refr_transcripts = agn_gene_locus_num_refr_transcripts(locus);
            locus_summaries[j].pred_transcripts = agn_gene_locus_num_pred_transcripts(locus);
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
                       rank, agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
                       agn_gene_locus_get_end(locus), agn_gene_locus_get_length(locus),
                       agn_gene_locus_num_refr_transcripts(locus),
                       agn_gene_locus_num_pred_transcripts(locus) );
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

            PeCompEvaluation data;
            //agn_gene_locus_get_summary_data(locus, &data);
            agn_gene_locus_aggregate_results(locus, &data);
            agn_gene_locus_aggregate_results(locus, &summary_data_local);
            GtArray *reported_pairs = agn_gene_locus_find_best_pairs(locus);
            if(reported_pairs == NULL)
              locus_summaries[j].reported = 0;
            else
              locus_summaries[j].reported = gt_array_size(reported_pairs);
            locus_summaries[j].counts = data.counts;
          }

          if(!options.summary_only)
          {
            //PeCompEvaluation data;
            //agn_gene_locus_get_summary_data(locus, &data);
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
              AgnGeneLocusPngMetadata metadata;
              pe_gene_locus_get_png_filename(locus, metadata.filename, options.outfilename);
              metadata.graphic_width = pe_gene_locus_get_graphic_width(locus);
              sprintf(metadata.stylefile, "%s/pe.style", options.data_path);
              metadata.refrfile = options.refrfile;
              metadata.predfile = options.predfile;
              metadata.refrlabel = options.refrlabel;
              metadata.predlabel = options.predlabel;
              metadata.track_order_func = pe_track_order;

              agn_gene_locus_print_png(locus, &metadata);
            }
#endif
          }

//          if(options.locusgff3)
//            agn_gene_locus_to_gff3(locus, locusgff3file);
        }
// fprintf(stderr, "DELETEME two\n");
        agn_gene_locus_delete(locus);
// fprintf(stderr, "DELETEME three\n");
      } // End parallel for loop

      #pragma omp critical(aggregate_stats)
      {
        pe_comp_evalutation_combine(&summary_data, &summary_data_local);
        pe_comp_evalutation_combine(&seqlevel_summary_data[i], &summary_data_local);
      }
    } // End pragma omp parallel

    if(!options.summary_only)
    {
      if(options.html)
      {
        for(j = 0; j < gt_array_size(seq_loci); j++)
        {
          AgnGeneLocusSummary *locus_summary = locus_summaries + j;
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
  agn_logger_delete(logger);
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
