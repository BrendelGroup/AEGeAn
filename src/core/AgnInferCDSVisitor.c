#include "core/array_api.h"
//#include "extended/node_visitor_api.h"
#include "AgnInferCDSVisitor.h"

#define infer_cds_visitor_cast(GV)\
        gt_node_visitor_cast(infer_cds_visitor_class(), GV)

//------------------------------------------------------------------------------
// Data structure definitions
//------------------------------------------------------------------------------

struct AgnInferCDSVisitor
{
  const GtNodeVisitor parent_instance;
  GtArray *cds;
  GtArray *utrs;
  GtArray *exons;
  GtArray *starts;
  GtArray *stops;
  GtUword cdscounter;
  GtLogger *logger;
};


//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function If the mRNA's CDS is discontinuous, ensure each CDS feature is
 * labeled as a multifeature.
 */
static void infer_cds_visitor_check_cds_multi(AgnInferCDSVisitor *v);

/**
 * @function Check and correct phase attributes for each CDS feature.
 */
static void infer_cds_visitor_check_cds_phase(AgnInferCDSVisitor *v);

/**
 * @function If start codon is provided explicitly, ensure it agrees with CDS,
 * whether the CDS is provided explicitly or implicitly inferred. If start codon
 * is not provided explicitly, infer it from CDS if possible.
 */
static void infer_cds_visitor_check_start(AgnInferCDSVisitor *v);

/**
 * @function If stop codon is provided explicitly, ensure it agrees with CDS,
 * whether the CDS is provided explicitly or implicitly inferred. If stop codon
 * is not provided explicitly, infer it from CDS if possible.
 */
static void infer_cds_visitor_check_stop(AgnInferCDSVisitor *v);

/**
 * @function Implement the interface to the GtNodeVisitor class.
 */
static const GtNodeVisitorClass *infer_cds_visitor_class();

/**
 * @function Infer CDS for any mRNAs that have none specified but have exons and
 * start/stop codons explicitly specified.
 */
static void infer_cds_visitor_infer_cds(AgnInferCDSVisitor *v);

/**
 * @function Infer UTRs for any mRNAs that have none specified but do have exons
 * and start/stop codons and/or CDS explicitly specified.
 */
static void infer_cds_visitor_infer_utrs(AgnInferCDSVisitor *v);

/**
 * @function Identify any mRNA subfeatures associated with this top-level
 * feature and apply the CDS inference procedure.
 */
static int infer_cds_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn,
                                                GtError *error);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

GtNodeStream* agn_infer_cds_stream_new(GtNodeStream *in, GtLogger *logger)
{
  GtNodeVisitor *nv = agn_infer_cds_visitor_new(logger);
  GtNodeStream *ns = gt_visitor_stream_new(in, nv);
  return ns;
}

GtNodeVisitor *agn_infer_cds_visitor_new(GtLogger *logger)
{
  GtNodeVisitor *nv;
  nv = gt_node_visitor_create(infer_cds_visitor_class());
  AgnInferCDSVisitor *v = infer_cds_visitor_cast(nv);
  v->logger = logger;
  v->cdscounter = 0;
infer_cds_visitor_check_cds_multi(v);
infer_cds_visitor_check_cds_phase(v);
infer_cds_visitor_check_start(v);
infer_cds_visitor_check_stop(v);
infer_cds_visitor_infer_cds(v);
infer_cds_visitor_infer_utrs(v);
  return nv;
}

bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)
{
  return false;
}

static void infer_cds_visitor_check_cds_multi(AgnInferCDSVisitor *v)
{
  
}

static void infer_cds_visitor_check_cds_phase(AgnInferCDSVisitor *v)
{
  
}

static void infer_cds_visitor_check_start(AgnInferCDSVisitor *v)
{
  
}

static void infer_cds_visitor_check_stop(AgnInferCDSVisitor *v)
{
  
}

static const GtNodeVisitorClass *infer_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if(!nvc)
  {
    nvc = gt_node_visitor_class_new(sizeof (AgnInferCDSVisitor), NULL, NULL,
                                    infer_cds_visitor_visit_feature_node, NULL,
                                    NULL, NULL);
  }
  return nvc;
}

static void infer_cds_visitor_infer_cds(AgnInferCDSVisitor *v)
{
  
}

static void infer_cds_visitor_infer_utrs(AgnInferCDSVisitor *v)
{
  
}

static int infer_cds_visitor_visit_feature_node(GtNodeVisitor *nv,
                                                GtFeatureNode *fn,
                                                GtError *error)
{
  return 0;
}
