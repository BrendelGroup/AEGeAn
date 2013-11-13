#include "PeProcedure.h"
#include "PeReports.h"

//------------------------------------------------------------------------------
// Prototype(s) for private function(s)
//------------------------------------------------------------------------------

/**
 * @function FIXME
 */
static void pe_check_filehandle_risk(GtUword numseqids);


//------------------------------------------------------------------------------
// Method/function implementations
//------------------------------------------------------------------------------

void pe_aggregate_results(AgnCompEvaluation *overall_eval,
                          GtArray **seqlevel_evalsp, GtArray *loci,
                          GtArray *seqfiles, GtHashmap *comp_evals,
                          GtHashmap *locus_summaries, PeOptions *options)
{
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin aggregating locus-level results\n", stderr);

  agn_comp_evaluation_init(overall_eval);
  GtArray *seqlevel_evals = gt_array_new( sizeof(AgnCompEvaluation) );
  GtUword i;
  for(i = 0; i < gt_array_size(seqfiles); i++)
  {
    FILE *seqfile = *(FILE **)gt_array_get(seqfiles, i);
    AgnCompEvaluation seqeval;
    agn_comp_evaluation_init(&seqeval);

    int j;
    GtArray *seqloci = *(GtArray **)gt_array_get(loci, i);
    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, j);
      AgnCompEvaluation *eval = gt_hashmap_get(comp_evals, locus);
      agn_comp_evaluation_combine(&seqeval, eval);
      agn_comp_evaluation_combine(overall_eval, eval);
      AgnGeneLocusSummary *locsum = gt_hashmap_get(locus_summaries, locus);
      // FIXME Should this block be placed elsewhere?
      if(strcmp(options->outfmt, "html") == 0 && !options->summary_only)
      {
        pe_print_locus_to_seqfile(seqfile, locsum->start, locsum->end,
                                  locsum->length, locsum->refrtrans,
                                  locsum->predtrans, &locsum->counts);
      }
    }
    gt_array_add(seqlevel_evals, seqeval);
    gt_array_delete(seqloci);
  }

  *seqlevel_evalsp = seqlevel_evals;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished aggregating locus-"
                          "level results (%ld.%06ld seconds)\n", stderr);
  gt_timer_delete(timer);
}

static void pe_check_filehandle_risk(GtUword numseqids)
{
  char *cmd = "ulimit -a | perl -ne 'if(m/\\(-n\\) (\\d+)/){ print $1 }'";
  char buffer[1024];
  FILE *stream = popen(cmd, "r");
  fgets(buffer, 1024, stream);
  fclose(stream);
  int openfilesallowed = atoi(buffer);
  if(openfilesallowed < numseqids + 1)
  {
    fprintf(stderr, "warning: max number of open files is %d, but there are "
            "%lu sequences to be compared; if ParsEval crashes, this is "
            "probably why; use 'ulimit -S -n $newlimit' to adjust this "
            "setting\n", openfilesallowed, numseqids);
  }
}

void pe_comparative_analysis(AgnLocusIndex *locusindex, GtHashmap **comp_evalsp,
                             GtHashmap **locus_summariesp, GtStrArray *seqids,
                             GtArray *seqfiles, GtArray *loci,
                             PeOptions *options)
{
  GtHashmap *comp_evals;
  GtHashmap *locus_summaries;

  AgnLogger *logger = agn_logger_new();
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin comparative analysis\n", stderr);
  locus_summaries = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);
  comp_evals = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree)gt_free_mem);

  int i, j;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    GtArray *seqloci = *(GtArray **)gt_array_get(loci, i);
    PeAnalysisData analysis_data;
    if(!options->summary_only)
      analysis_data.seqfile = *(FILE **)gt_array_get(seqfiles, i);
    analysis_data.options = options;

    for(j = 0; j < gt_array_size(seqloci); j++)
    {
      AgnGeneLocus *locus = *(AgnGeneLocus **)gt_array_get(seqloci, j);
      AgnCompEvaluation *compeval = gt_malloc( sizeof(AgnCompEvaluation) );
      agn_comp_evaluation_init(compeval);
      gt_hashmap_add(comp_evals, locus, compeval);
      AgnGeneLocusSummary *locsum = gt_malloc( sizeof(AgnGeneLocusSummary) );
      agn_gene_locus_summary_init(locsum);
      gt_hashmap_add(locus_summaries, locus, locsum);
    }
    analysis_data.comp_evals = comp_evals;
    analysis_data.locus_summaries = locus_summaries;

    agn_locus_index_comparative_analysis(locusindex, seqid,
                                     (AgnLocusIndexVisitFunc)pe_pre_analysis,
                                     (AgnLocusIndexVisitFunc)pe_post_analysis,
                                     &analysis_data, logger);
    if(options->debug)
      agn_logger_print_status(logger, stderr, "comparative analysis");
    if(agn_logger_has_error(logger))
    {
      agn_logger_print_error(logger, stderr, "issues with comparative analysis "
                             "of sequence '%s'", seqid);
      exit(1);
    }
  }

  *comp_evalsp = comp_evals;
  *locus_summariesp = locus_summaries;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished comparative "
                          "analysis (%ld.%06ld seconds)\n", stderr);
  agn_logger_delete(logger);
  gt_timer_delete(timer);
}

GtUword pe_load_and_parse_loci(AgnLocusIndex **locusindexp, GtArray **locip,
                               GtStrArray **seqidsp, PeOptions *options,
                               AgnLogger *logger)
{
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin loading data and parsing loci\n", stderr);

  // Load loci into memory
  AgnLocusIndex *locusindex = agn_locus_index_new(false);
  GtUword total = agn_locus_index_parse_pairwise_disk(locusindex,
                            options->refrfile, options->predfile,
                            &options->filters, logger);

  // Collect IDs of all sequences annotated by input files
  GtStrArray *seqids = agn_locus_index_seqids(locusindex);
  pe_check_filehandle_risk(gt_str_array_size(seqids));
  GtArray *loci = gt_array_new( sizeof(GtArray *) );
  int i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    const char *seqid = gt_str_array_get(seqids, i);
    pe_seqid_check(seqid, logger);

    GtArray *seq_loci = agn_locus_index_get(locusindex, seqid);
    gt_array_sort(seq_loci,(GtCompare)agn_gene_locus_array_compare);
    gt_array_add(loci, seq_loci);
  }

  *locusindexp = locusindex;
  *locip = loci;
  *seqidsp = seqids;

  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished loading data and "
                          "parsing loci (%ld.%06ld seconds)\n", stderr);
  gt_timer_delete(timer);
  return total;
}

void pe_post_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  AgnCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
  GtArray *pairs = agn_gene_locus_pairs_to_report(locus);
  if(gt_array_size(pairs) > 0)
    agn_gene_locus_aggregate_results(locus, compeval);

  AgnGeneLocusSummary *locsum = gt_hashmap_get(data->locus_summaries, locus);
  locsum->start = agn_gene_locus_get_start(locus);
  locsum->end = agn_gene_locus_get_end(locus);
  locsum->length = agn_gene_locus_get_length(locus);
  locsum->refrtrans = agn_gene_locus_num_refr_transcripts(locus);
  locsum->predtrans = agn_gene_locus_num_pred_transcripts(locus);
  locsum->reported = gt_array_size(pairs);
  locsum->counts = compeval->counts;

  if(!data->options->summary_only)
  {
    pe_gene_locus_print_results(locus, data->seqfile, data->options);

#ifndef WITHOUT_CAIRO
    if(data->options->locus_graphics)
    {
      AgnGeneLocusPngMetadata metadata;
      pe_gene_locus_get_png_filename(locus, metadata.filename,
                                     data->options->outfilename);
      metadata.graphic_width = pe_gene_locus_get_graphic_width(locus);
      sprintf(metadata.stylefile, "%s/pe.style", data->options->data_path);
      metadata.refrfile = data->options->refrfile;
      metadata.predfile = data->options->predfile;
      metadata.refrlabel = data->options->refrlabel;
      metadata.predlabel = data->options->predlabel;
      metadata.track_order_func = pe_track_order;

      agn_gene_locus_print_png(locus, &metadata);
    }
#endif
  }

  agn_gene_locus_delete(locus);
}

void pe_pre_analysis(AgnGeneLocus *locus, PeAnalysisData *data)
{
  GtUword npairs = agn_gene_locus_num_clique_pairs(locus);
  if(data->options->complimit != 0 && npairs > data->options->complimit)
  {
    if(data->options->debug)
    {
      fprintf(stderr, "debug: locus %s[%lu, %lu] contains %lu transcript (or "
              "transcript clique) pairs, exceeds the limit of %d, moving on\n",
              agn_gene_locus_get_seqid(locus), agn_gene_locus_get_start(locus),
              agn_gene_locus_get_end(locus), npairs, data->options->complimit);
    }
  }
  else
  {
    AgnCompEvaluation *compeval = gt_hashmap_get(data->comp_evals, locus);
    compeval->counts.num_loci++;
    compeval->counts.refr_genes += agn_gene_locus_num_refr_genes(locus);
    compeval->counts.pred_genes += agn_gene_locus_num_pred_genes(locus);
    GtUword rt = agn_gene_locus_num_refr_transcripts(locus);
    GtUword pt = agn_gene_locus_num_pred_transcripts(locus);
    compeval->counts.refr_transcripts += rt;
    compeval->counts.pred_transcripts += pt;
    if(agn_gene_locus_num_refr_genes(locus) == 0)
      compeval->counts.unique_pred++;
    else if(agn_gene_locus_num_pred_genes(locus) == 0)
      compeval->counts.unique_refr++;
  }
}

GtArray *pe_prep_output(GtStrArray *seqids, PeOptions *options)
{
  GtArray *seqfiles = gt_array_new( sizeof(FILE *) );

  if(strcmp(options->outfmt, "csv") == 0)
    pe_print_csv_header(options->outfile);

  GtUword i;
  for(i = 0; i < gt_str_array_size(seqids); i++)
  {
    FILE *seqfile = NULL;
    if(options->summary_only)
    {
      gt_array_add(seqfiles, seqfile);
    }
    else
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char filename[512];
      if(strcmp(options->outfmt, "html") == 0)
      {
        char dircmd[512];
        sprintf(dircmd, "mkdir %s/%s", options->outfilename, seqid);
        if(system(dircmd) != 0)
        {
          fprintf(stderr, "error: could not open directory '%s/%s'\n",
                  options->outfilename, seqid);
          exit(1);
        }
        if(options->debug)
          fprintf(stderr, "debug: opening directory '%s'\n", dircmd);

        sprintf(dircmd, "ln -s ../LICENSE %s/%s/LICENSE", options->outfilename,
                seqid);
        if(system(dircmd) != 0)
        {
          fputs("warning: could not create symbolic link to LICENSE\n", stderr);
        }

        // Create summary page for this sequence
        sprintf(filename, "%s/%s/index.html", options->outfilename, seqid);
        if(options->debug)
          fprintf(stderr, "debug: opening outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w", stderr);
        pe_print_seqfile_header(seqfile, seqid);
        gt_array_add(seqfiles, seqfile);
      }
      else
      {
        sprintf(filename, "%s.%s", options->outfilename, seqid);
        if(options->debug)
          fprintf(stderr, "debug: opening temp outfile '%s'\n", filename);
        seqfile = agn_fopen(filename, "w", stderr);
        gt_array_add(seqfiles, seqfile);
      }
    }
  }

  return seqfiles;
}

void pe_print_combine_output(GtStrArray *seqids, GtArray *seqfiles,
                             PeOptions *options)
{
  GtTimer *timer = gt_timer_new();
  gt_timer_start(timer);
  fputs("[ParsEval] Begin printing summary, combining output\n", stderr);

  if(!options->summary_only)
  {
    GtUword i;
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      FILE *seqfile = *(FILE **)gt_array_get(seqfiles, i);
      if(strcmp(options->outfmt, "html") == 0)
        pe_print_seqfile_footer(seqfile);
      fclose(seqfile);
    }
  }

  if(options->outfile == stdout)
    fflush(stdout);
  else
    fclose(options->outfile);
  if(!strcmp(options->outfmt, "html") == 0 && !options->summary_only)
  {
    int result;
    GtUword i;
    for(i = 0; i < gt_str_array_size(seqids); i++)
    {
      const char *seqid = gt_str_array_get(seqids, i);
      char command[1024];
      sprintf(command, "cat %s.%s >> %s && rm %s.%s", options->outfilename,
              seqid, options->outfilename, options->outfilename, seqid );
      if(options->outfile == stdout)
      {
        sprintf(command, "cat %s.%s && rm %s.%s", options->outfilename, seqid,
                options->outfilename, seqid );
      }
      if(options->debug)
        fprintf(stderr, "debug: merging output files: %s\n", command);
      result = system(command);
      if(result)
      {
        fprintf(stderr, "[ParsEval] error: issue merging GFF3 files: %s\n",
                command);
        exit(1);
      }
    }
  }
  gt_timer_stop(timer);
  gt_timer_show_formatted(timer, "[ParsEval] Finished printing summary, "
                          "combining output (%ld.%06ld seconds)\n", stderr);
  if(options->outfile == stdout)
    fclose(stdout);

  gt_timer_delete(timer);
}
