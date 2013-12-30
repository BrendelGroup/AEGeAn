AEGeAn C API
============

The AEGeAn Toolkit relies heavily on data
types implemented by the GenomeTools library. For data types beginning with
``Gt``, see the GenomeTools API documentation at
http://genometools.org/libgenometools.html.

Class AgnCliquePair
-------------------

.. c:type:: AgnCliquePair

  The AgnCliquePair class facilitates comparison of two alternative sources of annotation for the same sequence. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCliquePair.h>`_.

.. c:function:: AgnCompClassification agn_clique_pair_classify(AgnCliquePair *pair)

  Based on the already-computed comparison statistics, classify this clique pair as a perfect match, a CDS match, etc. See :c:type:`AgnCompClassification`.

.. c:function:: void agn_clique_pair_comparison_aggregate(AgnCliquePair *pair, AgnComparison *comp)

  Add this clique pair's internal comparison stats to a larger set of aggregate stats.

.. c:function:: int agn_clique_pair_compare(void *p1, void *p2)

  Same as c:func:`agn_clique_pair_compare_direct`, but with pointer dereferencing.

.. c:function:: int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2)

  Determine which pair has higher comparison scores. Returns 1 if the first pair has better scores, -1 if the second pair has better scores, 0 if they are equal.

.. c:function:: int agn_clique_pair_compare_reverse(void *p1, void *p2)

  Negation of c:func:`agn_clique_pair_compare`.

.. c:function:: void agn_clique_pair_delete(AgnCliquePair *pair)

  Class destructor.

.. c:function:: AgnCliquePair* agn_clique_pair_new(AgnTranscriptClique *refr, AgnTranscriptClique *pred)

  Class constructor.

.. c:function:: bool agn_clique_pair_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Module AgnComparison
--------------------

Data structures and functions related to comparative assessment of gene/transcript annotations. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnComparison.h>`_.

.. c:type:: AgnCompStatsBinary

  This struct is used to aggregate counts and statistics regarding the structural-level comparison (i.e., at the level of whole CDS segments, whole exons, and whole UTRs) and analysis of gene structure. See header file for details.



.. c:type:: AgnCompStatsScaled

  This struct is used to aggregate counts and statistics regarding the nucleotide-level comparison and analysis of gene structure. See header file for details.



.. c:type:: AgnComparison

  This struct aggregates all the counts and stats that go into a comparison, including structural-level and nucleotide-level counts and stats. See header file for details.



.. c:type:: AgnCompClassification

  This enumerated type refers to all the possible outcomes when annotations from two different sources are compared: ``AGN_COMP_CLASS_UNCLASSIFIED``, ``AGN_COMP_CLASS_PERFECT_MATCH``, ``AGN_COMP_CLASS_MISLABELED``, ``AGN_COMP_CLASS_CDS_MATCH``, ``AGN_COMP_CLASS_EXON_MATCH``, ``AGN_COMP_CLASS_UTR_MATCH``, and ``AGN_COMP_CLASS_NON_MATCH``.



.. c:function:: void agn_comparison_aggregate(AgnComparison *agg_cmp, AgnComparison *cmp)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comparison_init(AgnComparison *comparison)

  Initialize comparison stats to default values.

.. c:function:: void agn_comparison_resolve(AgnComparison *comparison)

  Calculate stats from the given counts.

.. c:function:: void agn_comp_stats_binary_aggregate(AgnCompStatsBinary *agg_stats, AgnCompStatsBinary *stats)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comp_stats_binary_init(AgnCompStatsBinary *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats)

  Calculate stats from the given counts.

.. c:function:: void agn_comp_stats_scaled_aggregate(AgnCompStatsScaled *agg_stats, AgnCompStatsScaled *stats)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats)

  Calculate stats from the given counts.

Class AgnTranscriptClique
-------------------------

.. c:type:: AgnTranscriptClique

  The purpose of the AgnTranscriptClique class is to store data pertaining to an individual maximal transcript clique. This clique may only contain a single transcript, or it may contain many. The only stipulation is that the transcripts do not overlap. Under the hood, each ``AgnTranscriptClique`` instance is a pseudo node (a GtFeatureNode object) with one or more transcript features as direct children. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTranscriptClique.h>`_.

.. c:type:: typedef void (*AgnCliqueVisitFunc)(GtFeatureNode*, void*)

   The signature that functions must match to be applied to each transcript in the given clique. The function will be called once for each transcript in the clique. The transcript will be passed as the first argument, and a second argument is available for an optional pointer to supplementary data (if needed). See :c:func:`agn_transcript_clique_traverse`.

.. c:function:: void agn_transcript_clique_add(AgnTranscriptClique *clique, GtFeatureNode *transcript)

  Add a transcript to this clique.

.. c:function:: GtUword agn_transcript_clique_cds_length(AgnTranscriptClique *clique)

  Get the combined CDS length (in base pairs) for all transcripts in this clique.

.. c:function:: AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)

  Make a shallow copy of this transcript clique.

.. c:function:: void agn_transcript_clique_delete(AgnTranscriptClique *clique)

  Class destructor.

.. c:function:: const char *agn_transcript_clique_get_model_vector(AgnTranscriptClique *clique)

  Get a pointer to the string representing this clique's transcript structure.

.. c:function:: bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique, GtHashmap *map)

  Determine whether any of the transcript IDs associated with this clique are keys in the given hash map.

.. c:function:: const char *agn_transcript_clique_id(AgnTranscriptClique *clique)

  Retrieve the ID attribute of the transcript associated with this clique.

.. c:function:: AgnTranscriptClique *agn_transcript_clique_new(AgnSequenceRegion *region)

  Class constructor. ``locusrange`` should be a pointer to the genomic coordinates of the locus to which this transcript clique belongs.

.. c:function:: GtUword agn_transcript_clique_num_exons(AgnTranscriptClique *clique)

  Get the number of exons in this clique.

.. c:function:: GtUword agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)

  Get the number of UTR segments in this clique.

.. c:function:: void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique, GtHashmap *map)

  Add all of the IDs associated with this clique to the given hash map.

.. c:function:: GtUword agn_transcript_clique_size(AgnTranscriptClique *clique)

  Get the number of transcripts in this clique.

.. c:function:: GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique)

  Get an array containing all the transcripts in this clique. User is responsible for deleting the array.

.. c:function:: void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream, const char *prefix)

  Print the transcript clique to the given outstream in GFF3 format, optionally with a prefix.

.. c:function:: void agn_transcript_clique_traverse(AgnTranscriptClique *clique, AgnCliqueVisitFunc func, void *funcdata)

  Apply ``func`` to each transcript in the clique. See :c:type:`AgnCliqueVisitFunc`.

.. c:function:: bool agn_transcript_clique_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Module AgnTypecheck
-------------------

Functions for testing feature types. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTypecheck.h>`_.

.. c:function:: bool agn_typecheck_cds(GtFeatureNode *fn)

  Returns true if the given feature is a CDS; false otherwise.

.. c:function:: bool agn_typecheck_exon(GtFeatureNode *fn)

  Returns true if the given feature is an exon; false otherwise.

.. c:function:: bool agn_typecheck_intron(GtFeatureNode *fn)

  Returns true if the given feature is an intron; false otherwise.

.. c:function:: bool agn_typecheck_mrna(GtFeatureNode *fn)

  Returns true if the given feature is an mRNA; false otherwise.

.. c:function:: bool agn_typecheck_transcript(GtFeatureNode *fn)

  Returns true if the given feature is an mRNA, tRNA, or rRNA; false otherwise.

.. c:function:: bool agn_typecheck_utr(GtFeatureNode *fn)

  Returns true if the given feature is a UTR; false otherwise.

.. c:function:: bool agn_typecheck_utr3p(GtFeatureNode *fn)

  Returns true if the given feature is a 3' UTR; false otherwise.

.. c:function:: bool agn_typecheck_utr5p(GtFeatureNode *fn)

  Returns true if the given feature is a 5' UTR; false otherwise.

Class AgnUnitTest
-----------------

.. c:type:: AgnUnitTest

  Class used for unit testing of classes and modules. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUnitTest.h>`_.

.. c:function:: void agn_unit_test_delete(AgnUnitTest *test)

  Destructor.

.. c:function:: AgnUnitTest *agn_unit_test_new(const char *label, bool (*testfunc)(AgnUnitTest *))

  Class constructor, where ``label`` is a label for the test and ``testfunc`` is a pointer to the function that will execute the test.

.. c:function:: void agn_unit_test_print(AgnUnitTest *test, FILE *outstream)

  Prints results of the unit test to ``outstream``.

.. c:function:: void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success)

  Add a result to this unit test.

.. c:function:: void agn_unit_test_run(AgnUnitTest *test)

  Run the unit test.

