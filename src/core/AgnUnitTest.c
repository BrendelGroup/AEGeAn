#include "AgnUnitTest.h"
#include "AgnUtils.h"

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnUnitTest
{
  char *label;
  bool (*testfunc)(struct AgnUnitTest *);
  GtArray *results;
  bool passed;
};

typedef struct
{
  char *label;
  bool success;
} UnitTestResult;


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_unit_test_delete(AgnUnitTest *test)
{
  gt_free(test->label);
  while(gt_array_size(test->results) > 0)
  {
    UnitTestResult *result = gt_array_pop(test->results);
    gt_free(result->label);
  }
  gt_array_delete(test->results);
  gt_free(test);
}

AgnUnitTest *agn_unit_test_new(const char *label,
                               bool (*testfunc)(AgnUnitTest *))
{
  AgnUnitTest *test = gt_malloc( sizeof(AgnUnitTest) );
  test->label = gt_cstr_dup(label);
  test->testfunc = testfunc;
  test->results = gt_array_new( sizeof(UnitTestResult) );
  test->passed = false;
  return test;
}

void agn_unit_test_print(AgnUnitTest *test, FILE *outstream)
{
  const char *successstr = "FAILURE";
  if(test->passed)
    successstr = "SUCCESS";
  fprintf(outstream, "    %-42s | %s\n", test->label, successstr);

  GtUword i;
  for(i = 0; i < gt_array_size(test->results); i++)
  {
    UnitTestResult *result = gt_array_get(test->results, i);
    const char *resultstr = "FAIL";
    if(result->success)
      resultstr = "PASS";
    fprintf(outstream, "        | %-36s | %s\n", result->label, resultstr);
  }
}

void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success)
{
  UnitTestResult result;
  result.label = gt_cstr_dup(label);
  result.success = success;
  gt_array_add(test->results, result);
}

bool agn_unit_test_success(AgnUnitTest *test)
{
  agn_assert(gt_array_size(test->results) > 0);
  GtUword i;
  for(i = 0; i < gt_array_size(test->results); i++)
  {
    UnitTestResult *result = gt_array_get(test->results, i);
    if(!result->success)
      return false;
  }
  return true;
}

void agn_unit_test_run(AgnUnitTest *test)
{
  test->passed = test->testfunc(test);
}
