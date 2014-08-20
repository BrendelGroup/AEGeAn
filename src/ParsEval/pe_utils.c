#include <stdio.h>
#include <string.h>
#include <time.h>
#include "pe_options.h"
#include "pe_utils.h"

char *pe_get_start_time()
{
  time_t start_time;
  struct tm *start_time_info;
  time(&start_time);
  start_time_info = localtime(&start_time);

  char timestr[128];
  strftime(timestr, 128, "%d %b %Y, %I:%M%p", start_time_info);
  return gt_cstr_dup(timestr);
}

void pe_summary_html_overview(FILE *outstream, void *data)
{
  int x;
  PeHtmlOverviewData *odata = data;
  fprintf(outstream,
          "      <pre class=\"command\">\n"
          "Started:                %s\n"
          "Reference annotations:  %s\n"
          "Prediction annotations: %s\n"
          "Executing command:      ",
          odata->start_time, odata->refrlabel, odata->predlabel);
  for(x = 0; x < odata->argc; x++)
  {
    fprintf(outstream, "%s ", odata->argv[x]);
  }
  fprintf(outstream, "</pre>\n\n");
}

void pe_summary_header(ParsEvalOptions *options, FILE *outstream,
                       char *start_time, int argc, char **argv)
{
  fprintf(outstream,
          "============================================================\n"
          "========== ParsEval Summary\n"
          "============================================================\n\n");

  fprintf(outstream, "Started:                %s\n", start_time);

  if(options->refrlabel != NULL)
    fprintf(outstream, "Reference annotations:  %s\n", options->refrlabel);
  else
    fprintf(outstream, "Reference annotations:  %s\n", options->refrfile);
  if(options->predlabel != NULL)
    fprintf(outstream, "Prediction annotations: %s\n", options->predlabel);
  else
    fprintf(outstream, "Prediction annotations: %s\n", options->predfile);
  fprintf(outstream, "Executing command:      ");

  int x;
  for(x = 0; x < argc; x++)
  {
    fprintf(outstream, "%s ", argv[x]);
  }
  fprintf(outstream, "\n\n");
}
