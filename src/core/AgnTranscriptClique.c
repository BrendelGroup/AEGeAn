#include "AgnTranscriptClique.h"
#include "AgnTypecheck.h"

//------------------------------------------------------------------------------
// Prototypes for private functions
//------------------------------------------------------------------------------

/**
 * @function Traversal function for determining the combined length of CDS(s)
 * for the mRNA(s) associated with this clique. Unit of length is amino acids.
 */
static void clique_cds_length(GtFeatureNode *fn, GtWord *length);

/**
 * @function Traversal function for copying the contents of one clique to
 * another.
 */
static void clique_copy(GtFeatureNode *fn, AgnTranscriptClique *newclique);

/**
 * @function Traversal function for determining the number of exons associated
 * with this transcript clique.
 */
static void clique_exon_count(GtFeatureNode *fn, GtWord *count);

/**
 * @function Traversal function for putting IDs of each transcript in this
 * clique in the given hash map.
 */
static void clique_ids_put(GtFeatureNode *fn, GtHashmap *map);

/**
 * @function Traversal function for determining the number of transcripts
 * associated with this transcript clique.
 */
static void clique_size(GtFeatureNode *fn, GtWord *count);

/**
 * @function Traversal function for copying members of this clique to an array.
 */
static void clique_to_array(GtFeatureNode *fn, GtArray *transcripts);

/**
 * @function Traversal function for printing members of this clique to an output
 * file in GFF3 format.
 */
static void clique_to_gff3(GtFeatureNode *fn, GtNodeVisitor *nv);

/**
 * @function Traverse every feature in this transcript clique's feature graph,
 * including all descendents of all transcripts, and apply ``func`` to each
 * feature (with ``funcdata`` available for a pointer to auxiliary data).
 */
static void clique_traverse(AgnTranscriptClique *clique,
                            AgnCliqueVisitFunc func,
                            void *funcdata);

/**
 * @function Traverse only the direct children of the clique pseudo-node (the
 * transcripts themselves) and apply ``func`` to each transcript (with
 * ``funcdata`` available for a pointer to auxiliary data).
 */
static void clique_traverse_direct(AgnTranscriptClique *clique,
                                   AgnCliqueVisitFunc func,
                                   void *funcdata);

/**
 * @function Traversal function for determining the number of UTR segments
 * associated with this transcript clique.
 */
static void clique_utr_count(GtFeatureNode *fn, GtWord *count);

/**
 * @function Update the clique's model vector whenever a new transcript is
 * added.
 */
static void clique_vector_update(AgnTranscriptClique *clique,
                                 GtFeatureNode *transcript);


//------------------------------------------------------------------------------
// Method implementations
//------------------------------------------------------------------------------

void agn_transcript_clique_add(AgnTranscriptClique *clique,
                               GtFeatureNode *transcript)
{
  gt_assert(agn_typecheck_transcript(transcript));
  GtFeatureNode *cliquefn = gt_feature_node_cast(clique);
  gt_feature_node_add_child(cliquefn, transcript);
  clique_vector_update(clique, transcript);
}

GtUword agn_transcript_clique_cds_length(AgnTranscriptClique *clique)
{
  GtUword length = 0;
  clique_traverse(clique, (AgnCliqueVisitFunc)clique_cds_length, &length);
  if(length % 3 != 0)
  {
    // FIXME
    // fprintf(stderr, "error: CDS length not divisible by 3\n");
  }
  return length / 3;
}

AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)
{
  GtStr *seqid = gt_genome_node_get_seqid(clique);
  GtRange range = gt_genome_node_get_range(clique);
  AgnSequenceRegion region = { gt_str_get(seqid), range };
  AgnTranscriptClique *newclique = agn_transcript_clique_new(&region);
  clique_traverse_direct(clique, (AgnCliqueVisitFunc)clique_copy, newclique);
  return newclique;
}

void agn_transcript_clique_delete(AgnTranscriptClique *clique)
{
  gt_genome_node_delete(clique);
}

const char *agn_transcript_clique_get_model_vector(AgnTranscriptClique *clique)
{
  return gt_genome_node_get_user_data(clique, "modelvector");
}

bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique,
                                          GtHashmap *map)
{
  bool idfound = false;

  GtFeatureNode *cliquefn = gt_feature_node_cast(clique);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(cliquefn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(current));
    const char *tid = gt_feature_node_get_attribute(current, "ID");
    if(gt_hashmap_get(map, tid) != NULL)
    {
      idfound = true;
      break;
    }
  }
  gt_feature_node_iterator_delete(iter);

  return idfound;
}

const char *agn_transcript_clique_id(AgnTranscriptClique *clique)
{
  char id[1024];
  char *idptr = id;
  unsigned count = 0;

  GtFeatureNode *cliquefn = gt_feature_node_cast(clique);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(cliquefn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    gt_assert(agn_typecheck_transcript(current));
    count++;
    if(count > 1)
      idptr += sprintf(idptr, ",");
    idptr += sprintf(idptr, "%s", gt_feature_node_get_attribute(current, "ID"));
  }
  gt_feature_node_iterator_delete(iter);

  return gt_cstr_dup(id);
}

AgnTranscriptClique *agn_transcript_clique_new(AgnSequenceRegion *region)
{
  GtStr *seqid = gt_str_new_cstr(region->seqid);
  AgnTranscriptClique *clique = gt_feature_node_new_pseudo(seqid,
                                                           region->range.start,
                                                           region->range.end,
                                                           GT_STRAND_BOTH);

  GtUword length = region->range.end - region->range.start + 1;
  char *modelvector = gt_malloc( sizeof(char) * (length + 1) );
  GtUword i;
  for(i = 0; i < length; i++)
    modelvector[i] = 'G';
  modelvector[length] = '\0';
  gt_genome_node_add_user_data(clique, "modelvector", modelvector,
                               gt_free_func);

  return clique;
}

GtUword agn_transcript_clique_num_exons(AgnTranscriptClique *clique)
{
  // FIXME Handle exons with same coordinates/strand
  GtUword count = 0;
  clique_traverse(clique, (AgnCliqueVisitFunc)clique_exon_count, &count);
  return count;
}

GtUword agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)
{
  // FIXME Handle UTR segments with same coordinates/strand
  GtUword count = 0;
  clique_traverse(clique, (AgnCliqueVisitFunc)clique_utr_count, &count);
  return count;
}

void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique,
                                           GtHashmap *map)
{
  clique_traverse_direct(clique, (AgnCliqueVisitFunc)clique_ids_put, map);
}

GtUword agn_transcript_clique_size(AgnTranscriptClique *clique)
{
  GtUword count = 0;
  clique_traverse_direct(clique, (AgnCliqueVisitFunc)clique_size, &count);
  return count;
}

GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique)
{
  GtArray *trans = gt_array_new( sizeof(GtFeatureNode *) );
  clique_traverse_direct(clique, (AgnCliqueVisitFunc)clique_to_array, trans);
  return trans;
}

void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream,
                                   const char *prefix)
{
  GtFile *outfile = gt_file_new_from_fileptr(outstream);
  GtNodeVisitor *nv = gt_gff3_visitor_new(outfile);
  // Never merged this patch with main GenomeTools repo
  // gt_gff3_visitor_set_output_prefix((GtGFF3Visitor *)nv, prefix);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor *)nv);
  clique_traverse_direct(clique, (AgnCliqueVisitFunc)clique_to_gff3, nv);
  gt_node_visitor_delete(nv);
  gt_file_delete_without_handle(outfile);
}

void agn_transcript_clique_traverse(AgnTranscriptClique *clique,
                                    AgnCliqueVisitFunc func,
                                    void *funcdata)
{
  clique_traverse_direct(clique, func, funcdata);
}

bool agn_transcript_clique_unit_test(AgnUnitTest *test)
{
  return false;
}

static void clique_cds_length(GtFeatureNode *fn, GtWord *length)
{
  if(agn_typecheck_cds(fn))
  {
    GtRange range = gt_genome_node_get_range((GtGenomeNode *)fn);
    (*length) += gt_range_length(&range);
  }
}

static void clique_copy(GtFeatureNode *fn, AgnTranscriptClique *newclique)
{
  agn_transcript_clique_add(newclique, fn);
}

static void clique_exon_count(GtFeatureNode *fn, GtWord *count)
{
  if(agn_typecheck_exon(fn))
    (*count)++;
}

static void clique_ids_put(GtFeatureNode *fn, GtHashmap *map)
{
  const char *tid = gt_feature_node_get_attribute(fn, "ID");
  gt_assert(gt_hashmap_get(map, tid) == NULL);
  gt_hashmap_add(map, (char *)tid, (char *)tid);
}

static void clique_size(GtFeatureNode *fn, GtWord *count)
{
  gt_assert(agn_typecheck_transcript(fn));
  (*count)++;
}

static void clique_to_array(GtFeatureNode *fn, GtArray *transcripts)
{
  gt_assert(agn_typecheck_transcript(fn));
  gt_array_add(transcripts, fn);
}

static void clique_to_gff3(GtFeatureNode *fn, GtNodeVisitor *nv)
{
  GtNodeVisitor *visitor = nv;
  GtError *error = gt_error_new();
  gt_genome_node_accept((GtGenomeNode *)fn, visitor, error);
  if(gt_error_is_set(error))
  {
    // FIXME?
    fprintf(stderr, "error with GFF3 output: %s\n", gt_error_get(error));
  }
  gt_error_delete(error);
}

static void clique_traverse(AgnTranscriptClique *clique,
                            AgnCliqueVisitFunc func,
                            void *funcdata)
{
  gt_assert(func);
  GtFeatureNode *cliquefn = gt_feature_node_cast(clique);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(cliquefn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    func(current, funcdata);
  }
  gt_feature_node_iterator_delete(iter);
}

static void clique_traverse_direct(AgnTranscriptClique *clique,
                                   AgnCliqueVisitFunc func,
                                   void *funcdata)
{
  gt_assert(func);
  GtFeatureNode *cliquefn = gt_feature_node_cast(clique);
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new_direct(cliquefn);
  GtFeatureNode *current;
  for(current = gt_feature_node_iterator_next(iter);
      current != NULL;
      current = gt_feature_node_iterator_next(iter))
  {
    func(current, funcdata);
  }
  gt_feature_node_iterator_delete(iter);
}

static void clique_utr_count(GtFeatureNode *fn, GtWord *count)
{
  if(agn_typecheck_utr(fn))
    (*count)++;
}

static void clique_vector_update(AgnTranscriptClique *clique,
                                 GtFeatureNode *transcript)
{
  GtRange locusrange = gt_genome_node_get_range(clique);
  char *modelvector = gt_genome_node_get_user_data(clique, "modelvector");
  GtFeatureNode *fn;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
  for(fn = gt_feature_node_iterator_next(iter);
      fn != NULL;
      fn = gt_feature_node_iterator_next(iter))
  {
    char c;
    if(agn_typecheck_cds(fn))
      c = 'C';
    else if(agn_typecheck_utr(fn))
    {
      gt_assert(agn_typecheck_utr3p(fn) || agn_typecheck_utr5p(fn));
      if(agn_typecheck_utr5p(fn))
        c = 'F';
      else
        c = 'T';
    }
    else if(agn_typecheck_intron(fn))
      c = 'I';
    else
      c = 'G';

    if(c != 'G')
    {
      GtUword fn_start = gt_genome_node_get_start((GtGenomeNode *)fn);
      GtUword fn_end = gt_genome_node_get_end((GtGenomeNode *)fn);
      GtUword i;

      for(i = fn_start - locusrange.start;
          i < fn_end - locusrange.start + 1;
          i++)
      {
        modelvector[i] = c;
      }
    }
  }
  gt_feature_node_iterator_delete(iter);
}
