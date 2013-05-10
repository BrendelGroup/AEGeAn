#include <string.h>
#include "genometools.h"
#include "VangDegreeConstraint.h"

DegreeConstraint vang_degree_constraint_parse(const char *degree_string)
{
  DegreeContext context;
  DegreeOperator op;
  unsigned int degree;

  char *degstrcpy = gt_cstr_dup(degree_string);
  char *contextstr = strtok(degstrcpy, " ");
  char *operatorstr = strtok(NULL, " ");
  char *degreestr = strtok(NULL, " ");

  if(strcmp(contextstr, "in") == 0)
    context = DEGREE_IN;
  else if(strcmp(contextstr, "out") == 0)
    context = DEGREE_OUT;
  else
  {
    fprintf( stderr, "error: unknown degree context '%s' (%s)\n", contextstr,
             degree_string );
    exit(1);
  }

  if(strcmp(operatorstr, "+") == 0)
    op = AT_LEAST;
  else if(strcmp(operatorstr, "-") == 0)
    op = AT_MOST;
  else if(strcmp(operatorstr, ">") == 0)
    op = GREATER_THAN;
  else if(strcmp(operatorstr, "<") == 0)
    op = LESS_THAN;
  else if(strcmp(operatorstr, "e") == 0)
    op = EQUALS;
  else
  {
    fprintf( stderr, "error: unknown degree constraint operator '%s' (%s)\n",
             operatorstr, degree_string );
    exit(1);
  }

  degree = atoi(degreestr);

  DegreeConstraint degree_constraint = {context, op, degree};
  gt_free(degstrcpy);

  return degree_constraint;
}

void vang_degree_constraint_to_string(DegreeConstraint *dc, FILE *outstream)
{
  const char *context = NULL;
  char operator = '\0';

  switch(dc->context)
  {
    case DEGREE_IN:
      context = "in";
      break;
    case DEGREE_OUT:
      context = "out";
      break;
    default:
      context = "undef";
      break;
  }

  switch(dc->operator)
  {
    case AT_LEAST:
      operator = '+';
      break;
    case AT_MOST:
      operator = '-';
      break;
    case GREATER_THAN:
      operator = '>';
      break;
    case LESS_THAN:
      operator = '<';
      break;
    case EQUALS:
      operator = 'e';
      break;
    default:
      operator = '?';
      break;
  }

  fprintf(outstream, "%s %c %u", context, operator, dc->degree);
}

