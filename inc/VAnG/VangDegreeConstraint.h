/**

Copyright (c) 2010-2014, Daniel S. Standage and CONTRIBUTORS

The AEGeAn Toolkit is distributed under the ISC License. See
the 'LICENSE' file in the AEGeAn source code distribution or
online at https://github.com/standage/AEGeAn/blob/master/LICENSE.

**/
#ifndef VANG_DEGREE_CONSTRAINT
#define VANG_DEGREE_CONSTRAINT

//----------------------------------------------------------------------------//
// Enumerated types and type definitions
//----------------------------------------------------------------------------//
enum DegreeContextEnum {DEGREE_IN, DEGREE_OUT};
enum DegreeOperatorEnum {AT_LEAST, AT_MOST, LESS_THAN, GREATER_THAN, EQUALS};
typedef struct DegreeConstraint DegreeConstraint;
typedef enum DegreeContextEnum DegreeContext;
typedef enum DegreeOperatorEnum DegreeOperator;


//----------------------------------------------------------------------------//
// Public data structure definition
//----------------------------------------------------------------------------//
struct DegreeConstraint
{
  DegreeContext context;
  DegreeOperator operator;
  unsigned int degree;
};


//----------------------------------------------------------------------------//
// Method prototypes
//----------------------------------------------------------------------------//

/**
 * Given a space-delimited string with 3 values, parse a degree constraint
 *
 * @param[in] degree_string    a space-delimited string with 3 values
 *                             corresponding to the context, operator, and
 *                             degree of a degree constraint
 * @returns                    a degree constraint object
 */
DegreeConstraint vang_degree_constraint_parse(const char *degree_string);

/**
 * Print a string representation of this degree constraint
 *
 * @param[in]  dc           degree constraint object
 * @param[out] outstream    output stream to which data will be printed
 */
void vang_degree_constraint_to_string(DegreeConstraint *dc, FILE *outstream);

#endif
