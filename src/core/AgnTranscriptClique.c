#include <string.h>
#include "AgnGeneLocus.h"
#include "AgnGtExtensions.h"
#include "AgnTestData.h"
#include "AgnTranscriptClique.h"
#include "AgnUnitTest.h"
#include "AgnUtils.h"

//----------------------------------------------------------------------------//
// Data structure definitions
//----------------------------------------------------------------------------//

/**
 * This data structure is pretty minimal. At one point it contained more than a
 * single member, and it may yet again in the future. For now though, there is
 * enough logic built around these data that a class is still necessary, even if
 * much of the core logic belongs to the GtDlist class.
 */
struct AgnTranscriptClique
{
  GtDlist *transcripts;
};

typedef struct
{
  GtHashmap *map;
  bool idfound;
} TranscriptIdCheckData;

typedef struct
{
  FILE *outstream;
  GtFeatureNode *first;
} PrintIdsData;


//----------------------------------------------------------------------------//
// Private method prototypes
//----------------------------------------------------------------------------//

/**
 * Traversal function for determining length of CDS(s) for the transcript(s) in
 * this clique.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] cdslength     pointer to an GtUword for storing the CDS length
 */
static void clique_cds_length(GtFeatureNode *transcript, void *cdslength);

/**
 * Traversal function for copying the contents of one clique to another.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] clique        new clique to which the transcript will be added
 */
static void clique_copy(GtFeatureNode *transcript, void *clique);

/**
 * Traversal function for checking whether the clique includes one or more
 * transcripts whose IDs are associated with the provided data.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] data          TranscriptIdCheckData, including a hashmap of IDs
 *                           and a boolean indicating whether the ID is in the
 *                           hashmap
 */
static void clique_id_check(GtFeatureNode *transcript, void *data);

/**
 * Traversal function for placing all transcript IDs associated with this clique
 * in the provided hashmap.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] map           hashmap of transcript IDs
 */
static void clique_ids_put(GtFeatureNode *transcript, void *map);

/**
 * Traversal function for determining the number of exons belonging to this
 * transcript clique.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] numexons      pointer to a GtUword for tracking exon #
 */
static void clique_num_exons(GtFeatureNode *transcript, void *numexons);

/**
 * Traversal function for determining the number of UTR segments belonging to
 * this transcript clique.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] numutrs       pointer to a GtUword for tracking UTR #
 */
static void clique_num_utrs(GtFeatureNode *transcript, void *numutrs);

/**
 * Print IDs to the given output stream.
 *
 * @param[in] transcript    transcript in the clique
 * @param[in] data          data necessary for printing the output
 */
static void clique_print_ids(GtFeatureNode *transcript, void *data);

/**
 * Traversal function for copying contents of this clique to an array.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] array         array to which transcripts will be added
 */
static void clique_to_array(GtFeatureNode *transcript, void *array);

/**
 * Traversal function for copying contents of this clique to an array.
 *
 * @param[in]  transcript    transcript in the clique
 * @param[out] nv            node visitor for generating GFF3 output
 */
static void clique_to_gff3(GtFeatureNode *transcript, void *nv);


//----------------------------------------------------------------------------//
// Method implementations
//----------------------------------------------------------------------------//

void agn_transcript_clique_add(AgnTranscriptClique *clique,
                               GtFeatureNode *transcript)
{
  gt_dlist_add(clique->transcripts, transcript);
}

GtUword agn_transcript_clique_cds_length(AgnTranscriptClique *clique)
{
  GtUword length = 0;
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_cds_length,
                                 &length);
  return length;
}

AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)
{
  AgnTranscriptClique *new = agn_transcript_clique_new();
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_copy, new);
  return new;
}

void agn_transcript_clique_delete(AgnTranscriptClique *clique)
{
  gt_dlist_delete(clique->transcripts);
  gt_free(clique);
  clique = NULL;
}

bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique,
                                          GtHashmap *map)
{
  TranscriptIdCheckData data = { map, false };
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_id_check,
                                 &data);
  return data.idfound;
}

const char *agn_transcript_clique_id(AgnTranscriptClique *clique)
{
  gt_assert(agn_transcript_clique_size(clique) == 1);
  GtDlistelem *elem = gt_dlist_first(clique->transcripts);
  GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
  return gt_feature_node_get_attribute(transcript, "ID");
}

AgnTranscriptClique* agn_transcript_clique_new()
{
  AgnTranscriptClique *clique = gt_malloc(sizeof(AgnTranscriptClique));
  clique->transcripts = gt_dlist_new((GtCompare)gt_genome_node_cmp);
  return clique;
}

GtUword agn_transcript_clique_num_exons(AgnTranscriptClique *clique)
{
  GtUword exon_count = 0;
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_num_exons,
                                 &exon_count);
  return exon_count;
}

GtUword agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)
{
  GtUword utr_count = 0;
  agn_transcript_clique_traverse(clique, clique_num_utrs, &utr_count);
  return utr_count;
}

void agn_transcript_clique_print_ids(AgnTranscriptClique *clique,
                                     FILE *outstream)
{
  if(gt_dlist_size(clique->transcripts) == 1)
  {
    GtDlistelem *elem = gt_dlist_first(clique->transcripts);
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    fprintf(outstream, "%s", gt_feature_node_get_attribute(transcript, "ID"));
    return;
  }
  fprintf(outstream, "[");
  GtDlistelem *first = gt_dlist_first(clique->transcripts);
  PrintIdsData data = { outstream, gt_dlistelem_get_data(first) };
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_print_ids,
                                 &data);
  fprintf(outstream, "]");

}

void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique,
                                           GtHashmap *map)
{
  agn_transcript_clique_traverse(clique, clique_ids_put, map);
}

GtUword agn_transcript_clique_size(AgnTranscriptClique *clique)
{
  return gt_dlist_size(clique->transcripts);
}

GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique)
{
  GtArray *new = gt_array_new( sizeof(GtFeatureNode *) );
  agn_transcript_clique_traverse(clique, (AgnCliqueVisitFunc)clique_to_array,
                                 new);
  return new;
}

void agn_transcript_clique_traverse(AgnTranscriptClique *clique,
                                    AgnCliqueVisitFunc func, void *funcdata)
{
  gt_assert(func);
  GtDlistelem *elem;
  for(elem = gt_dlist_first(clique->transcripts);
      elem != NULL;
      elem = gt_dlistelem_next(elem))
  {
    GtFeatureNode *transcript = gt_dlistelem_get_data(elem);
    func(transcript, funcdata);
  }
}

void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream,
                                   const char *prefix)
{
  GtFile *outfile = gt_file_new_from_fileptr(outstream);
  GtNodeVisitor *nv = gt_gff3_visitor_new(outfile);
  // Patch is not yet available from core GenomeTools repo
  // gt_gff3_visitor_set_output_prefix((GtGFF3Visitor *)nv, prefix);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor *)nv);
  agn_transcript_clique_traverse(clique,(AgnCliqueVisitFunc)clique_to_gff3, nv);
  gt_node_visitor_delete(nv);
  gt_file_delete_without_handle(outfile);
}

bool agn_transcript_clique_unit_test(AgnUnitTest *test)
{
  GtFeatureNode *eden = agn_test_data_eden();
  GtStr *seqid = gt_genome_node_get_seqid((GtGenomeNode *)eden);
  AgnGeneLocus *locus = agn_gene_locus_new(gt_str_get(seqid));
  agn_gene_locus_add_gene(locus, eden);

  GtArray *trans = agn_gene_locus_get_transcripts(locus);
  GtArray *cliques = agn_enumerate_feature_cliques(trans);
  bool parsearraypass = gt_array_size(cliques) == 3;
  agn_unit_test_result(test, "parse from array", parsearraypass);

  bool numtranspass = true;
  bool cdslenpass = true;
  bool iterpass = true;
  GtUword i;
  for(i = 0; i < gt_array_size(cliques); i++)
  {
    AgnTranscriptClique **tc = gt_array_get(cliques, i);
    if(agn_transcript_clique_size(*tc) != 1)
    {
      numtranspass = false;
    }
    const char *mrnaid = agn_transcript_clique_id(*tc);
    GtUword cdslength = agn_transcript_clique_cds_length(*tc);
    if(strcmp(mrnaid, "EDEN.1") == 0)
    {
      if(cdslength != 2305) cdslenpass = false;
    }
    else if(strcmp(mrnaid, "EDEN.2") == 0)
    {
      if(cdslength != 1402) cdslenpass = false;
    }
    else if(strcmp(mrnaid, "EDEN.3") == 0)
    {
      if(cdslength != 1704) cdslenpass = false;
    }
    else
    {
      iterpass = false;
    }
  }
  agn_unit_test_result(test, "transcripts per clique",      numtranspass);
  agn_unit_test_result(test, "clique CDS length",           cdslenpass);
  agn_unit_test_result(test, "iterate through transcripts", iterpass);

  agn_gene_locus_delete(locus);
  while(gt_array_size(cliques) > 0)
  {
    AgnTranscriptClique **tc = gt_array_pop(cliques);
    agn_transcript_clique_delete(*tc);
  }
  gt_array_delete(cliques);
  gt_array_delete(trans);
  gt_genome_node_delete((GtGenomeNode *)eden);

  return parsearraypass && numtranspass && cdslenpass && iterpass;
}

static void clique_cds_length(GtFeatureNode *transcript, void *cdslength)
{
  GtUword *length = cdslength;
  *length += agn_gt_feature_node_cds_length(transcript);
}

static void clique_copy(GtFeatureNode *transcript, void *clique)
{
  AgnTranscriptClique *cq = clique;
  agn_transcript_clique_add(cq, transcript);
}

static void clique_id_check(GtFeatureNode *transcript, void *data)
{
  TranscriptIdCheckData *dat = data;
  const char *tid = gt_feature_node_get_attribute(transcript, "ID");
  if(gt_hashmap_get(dat->map, tid) != NULL)
    dat->idfound = true;
}

static void clique_ids_put(GtFeatureNode *transcript, void *map)
{
  GtHashmap *mp = map;
  const char *tid = gt_feature_node_get_attribute(transcript, "ID");
  gt_assert(gt_hashmap_get(mp, tid) == NULL);
  gt_hashmap_add(map, (char *)tid, (char *)tid);
}

static void clique_num_exons(GtFeatureNode *transcript, void *numexons)
{
  GtUword *num = numexons;
  GtFeatureNode *current;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_gt_feature_node_is_exon_feature(current))
      *num += 1;
  }
  gt_feature_node_iterator_delete(iter);
}

static void clique_num_utrs(GtFeatureNode *transcript, void *numutrs)
{
  GtUword *num = numutrs;
  GtFeatureNode *current;
  GtFeatureNodeIterator *iter = gt_feature_node_iterator_new(transcript);
  for(current  = gt_feature_node_iterator_next(iter);
      current != NULL;
      current  = gt_feature_node_iterator_next(iter))
  {
    if(agn_gt_feature_node_is_utr_feature(current))
      *num += 1;
  }
  gt_feature_node_iterator_delete(iter);
}

static void clique_print_ids(GtFeatureNode *transcript, void *data)
{
  PrintIdsData *dat = data;
  const char *id = gt_feature_node_get_attribute(transcript, "ID");
  fprintf(dat->outstream, "%s", id);
  if(transcript != dat->first)
    fputc(',', dat->outstream);
}

static void clique_to_array(GtFeatureNode *transcript, void *array)
{
  GtArray *ar = (GtArray *)array;
  gt_array_add(ar, transcript);
}

static void clique_to_gff3(GtFeatureNode *transcript, void *nv)
{
  GtNodeVisitor *visitor = nv;
  GtError *error = gt_error_new();
  gt_genome_node_accept((GtGenomeNode *)transcript, visitor, error);
  if(gt_error_is_set(error))
    fprintf(stderr, "error with GFF3 output: %s\n", gt_error_get(error));
  gt_error_delete(error);
}
