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

.. c:function:: AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair)

  Return a pointer to the prediction annotation from this pair.

.. c:function:: AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair)

  Return a pointer to the reference annotation from this pair.

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

.. c:function:: void agn_comparison_print(AgnComparison *stats, FILE *outstream)

  Print the comparison stats to the given file.

.. c:function:: void agn_comparison_resolve(AgnComparison *comparison)

  Calculate stats from the given counts.

.. c:function:: bool agn_comparison_test(AgnComparison *c1, AgnComparison *c2)

  Returns true if c1 and c2 contain identical values, false otherwise.

.. c:function:: void agn_comp_stats_binary_aggregate(AgnCompStatsBinary *agg_stats, AgnCompStatsBinary *stats)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comp_stats_binary_init(AgnCompStatsBinary *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_binary_print(AgnCompStatsBinary *stats, FILE *outstream)

  Print the comparison stats to the given file.

.. c:function:: void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats)

  Calculate stats from the given counts.

.. c:function:: bool agn_comp_stats_binary_test(AgnCompStatsBinary *s1, AgnCompStatsBinary *s2)

  Returns true if s1 and s2 contain identical values, false otherwise.

.. c:function:: void agn_comp_stats_scaled_aggregate(AgnCompStatsScaled *agg_stats, AgnCompStatsScaled *stats)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_scaled_print(AgnCompStatsScaled *stats, FILE *outstream)

  Print the comparison stats to the given file.

.. c:function:: void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats)

  Calculate stats from the given counts.

.. c:function:: bool agn_comp_stats_scaled_test(AgnCompStatsScaled *s1, AgnCompStatsScaled *s2)

  Returns true if s1 and s2 contain identical values, false otherwise.

Class AgnFilterStream
---------------------

.. c:type:: AgnFilterStream

  Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream used to select features of a certain type from a node stream. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnFilterStream.h>`_.

.. c:function:: GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream, GtHashmap *typestokeep)

  Class constructor. The keys of the ``typestokeep`` hashmap should be the type(s) to be kept from the node stream. Any non-NULL value can be associated with those keys.

.. c:function:: bool agn_filter_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnInferCDSVisitor
------------------------

.. c:type:: AgnInferCDSVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring an mRNA's CDS from explicitly defined exon and start/stop codon features. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferCDSVisitor.h>`_.

.. c:function:: GtNodeStream* agn_infer_cds_stream_new(GtNodeStream *in, GtLogger *logger)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_infer_cds_visitor_new(GtLogger *logger)

  Constructor for the node visitor.

.. c:function:: bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnInferExonsVisitor
--------------------------

.. c:type:: AgnInferExonsVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring exon features when only CDS and UTR features are provided explicitly.  See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferExonsVisitor.h>`_.

.. c:function:: GtNodeStream* agn_infer_exons_stream_new(GtNodeStream *in, GtLogger *logger)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor* agn_infer_exons_visitor_new(GtLogger *logger)

  Class constructor for the node visitor.

.. c:function:: bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class.

Class AgnLocus
--------------

.. c:type:: AgnLocus

  The AgnLocus class represents gene loci and interval loci in memory and can be used to facilitate comparison of two different sources of annotation. Under the hood, each ``AgnLocus`` object is a feature node with one or more transcript features as direct children. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocus.h>`_.

.. c:type:: AgnComparisonSource

  When tracking the source of an annotation for comparison purposes, use this enumerated type to refer to reference (``REFERENCESOURCE``) vs prediction (``PREDICTIONSOURCE``) annotations. ``DEFAULTSOURCE`` is for when the source is not a concern.



.. c:type:: AgnLocusPngMetadata

  This data structure provides a convenient container for metadata needed to produce a PNG graphic for pairwise comparison loci.



.. c:function:: void agn_locus_add(AgnLocus *locus, GtFeatureNode *transcript, AgnComparisonSource source)

  Associate the given transcript annotation with this locus. Rather than calling this function directly, users are recommended to use one of the following macros: ``agn_locus_add_pred_transcript(locus, trans)`` and ``agn_locus_add_refr_transcript(locus, trans)``, to be used when keeping track of an annotation's source is important (i.e. for pairwise comparison); and ``agn_locus_add_transcript(locus, trans)`` otherwise.

.. c:function:: AgnLocus *agn_locus_clone(AgnLocus *locus)

  Do a semi-shallow copy of this data structure--for members whose data types support reference counting, the same pointer is used and the reference is incremented. For the other members a new object is created and populated with the same content.

.. c:function:: GtUword agn_locus_cds_length(AgnLocus *locus, AgnComparisonSource src)

  The combined length of all coding sequences associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_refr_cds_length(locus)`` for the combined length of all reference CDSs, ``agn_locus_pred_cds_length(locus)`` for the combined length of all prediction CDSs, and ``agn_locus_get_cds_length(locus)`` for the combined length of all CDSs.

.. c:function:: void agn_locus_comparative_analysis(AgnLocus *locus, GtUword maxtranscripts, GtUword maxpairs, GtLogger *logger)

  Compare every reference transcript clique with every prediction transcript clique. For gene loci with multiple transcript cliques, each comparison is not necessarily reported. Instead, we report the set of clique pairs that provides the optimal pairing of reference and prediction transcripts. If there are more reference transcript cliques than prediction cliques (or vice versa), these unmatched cliques are reported separately.

.. c:function:: int agn_locus_array_compare(const void *p1, const void *p2)

  Analog of ``strcmp`` for sorting AgnLocus objects. Loci are first sorted lexicographically by sequence ID, and then spatially by genomic coordinates.

.. c:function:: void agn_locus_comparison_aggregate(AgnLocus *locus, AgnComparison *comp)

  Add this locus' internal comparison stats to a larger set of aggregate stats.

.. c:function:: void agn_locus_delete(AgnLocus *locus)

  Class destructor.

.. c:function:: GtUword agn_locus_exon_num(AgnLocus *locus, AgnComparisonSource src)

  Get the number of exons for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_num_pred_exons(locus)`` for the number of prediction exons, ``agn_locus_num_refr_exons(locus)`` for the number of reference exons, or ``agn_locus_num_exons(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_locus_get_unique_pred_cliques(AgnLocus *locus)

  Get a list of all the prediction transcript cliques that have no corresponding reference transcript clique.

.. c:function:: GtArray *agn_locus_get_unique_refr_cliques(AgnLocus *locus)

  Get a list of all the reference transcript cliques that have no corresponding prediction transcript clique.

.. c:function:: AgnLocus* agn_locus_new(GtStr *seqid)

  Class constructor.

.. c:function:: GtUword agn_locus_num_clique_pairs(AgnLocus *locus)

  Return the number of clique pairs to be reported for this locus.

.. c:function:: void agn_locus_png_track_selector(GtBlock *block, GtStr *track,void *data)

  Track selector function for generating PNG graphics of pairwise comparison loci. The track name to will be written to ``track``.

.. c:function:: void agn_locus_print_png(AgnLocus *locus, AgnLocusPngMetadata *metadata)

  Print a PNG graphic for this locus.

.. c:function:: void agn_locus_print_transcript_mapping(AgnLocus *locus, FILE *outstream)

  Print a mapping of the transcript(s) associated with this locus in a two-column tab-delimited format: ``transcriptId<tab>locusId``.

.. c:function:: double agn_locus_splice_complexity(AgnLocus *locus, AgnComparisonSource src)

  Calculate the splice complexity of this gene locus. Rather than calling this method directly, users are recommended to use one of the following macros: ``agn_locus_prep_splice_complexity(locus)`` to calculate the splice complexity of just the prediction transcripts, ``agn_locus_refr_splice_complexity(locus)`` to calculate the splice complexity of just the reference transcripts, and ``agn_locus_calc_splice_complexity(locus)`` to calculate the splice complexity taking into account all transcripts.

.. c:function:: GtArray *agn_locus_transcripts(AgnLocus *locus, AgnComparisonSource src)

  Get the transcripts associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_transcripts(locus)`` to retrieve prediction transcripts, ``agn_locus_refr_transcripts(locus)`` to retrieve reference transcripts, or ``agn_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_locus_transcript_ids(AgnLocus *locus, AgnComparisonSource src)

  Get the transcript IDs associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_transcripts(locus)`` to retrieve prediction IDs, ``agn_locus_refr_transcripts(locus)`` to retrieve reference IDs, or ``agn_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_locus_transcript_num(AgnLocus *locus, AgnComparisonSource src)

  Get the number of transcripts for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_transcript_locus_num_pred_transcripts(locus)`` for the number of prediction transcripts, ``agn_transcript_locus_num_refr_transcripts(locus)`` for the number of reference transcripts, or ``agn_transcript_locus_num_transcripts(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: bool agn_locus_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

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

.. c:function:: bool agn_typecheck_gene(GtFeatureNode *fn)

  Returns true if the given feature is a gene; false otherwise.

.. c:function:: bool agn_typecheck_intron(GtFeatureNode *fn)

  Returns true if the given feature is an intron; false otherwise.

.. c:function:: bool agn_typecheck_mrna(GtFeatureNode *fn)

  Returns true if the given feature is an mRNA; false otherwise.

.. c:function:: GtArray *agn_typecheck_select(GtFeatureNode *fn, bool (*func)(GtFeatureNode *))

  Gather the children of a given feature that have a certain type. Type is tested by ``func``, which accepts a single ``GtFeatureNode`` object.

.. c:function:: bool agn_typecheck_start_codon(GtFeatureNode *fn)

  Returns true if the given feature is a start codon; false otherwise.

.. c:function:: bool agn_typecheck_stop_codon(GtFeatureNode *fn)

  Returns true if the given feature is a stop codon; false otherwise.

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

.. c:function:: bool agn_unit_test_success(AgnUnitTest *test)

  Returns true if all the results checked with this unit test passed, false otherwise.

.. c:function:: void agn_unit_test_run(AgnUnitTest *test)

  Run the unit test.

Module AgnUtils
---------------

Collection of assorted functions that are otherwise unrelated. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUtils.h>`_.

.. c:type:: AgnSequenceRegion

  This data structure combines sequence coordinates with a sequence ID to facilitate their usage together.



.. c:function:: GtArray* agn_array_copy(GtArray *source, size_t size)

  Similar to ``gt_array_copy``, except that array elements are treated as pointers and dereferenced before being added to the new array.

.. c:function:: double agn_calc_splice_complexity(GtArray *transcripts)

  Determine the splice complexity of the given set of transcripts.

.. c:function:: int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)

  Compare function for data type ``GtGenomeNode ``, needed for sorting ``GtGenomeNode `` stored in ``GtArray`` objects.

