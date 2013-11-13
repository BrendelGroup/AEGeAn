AEGeAn C API
============

The AEGeAn Toolkit relies heavily on data
types implemented by the GenomeTools library. For data types beginning with
``Gt``, see the GenomeTools API documentation at
http://genometools.org/libgenometools.html.

Class AgnCanonGeneStream
------------------------

.. c:type:: AgnCanonGeneStream

  Implements the GenomeTools ``GtNodeStream`` interfact. This is a node stream that returns only canonical protein-coding genes that pass stringent validation. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCanonGeneStream.h>`_.

.. c:function:: GtNodeStream* agn_canon_gene_stream_new(GtNodeStream *in_stream, AgnLogger *logger)

  Class constructor.

Class AgnCliquePair
-------------------

.. c:type:: AgnCliquePair

  The purpose of the AgnCliquePair class is to associate all of the data needed to compare transcripts from two alternative sources of annotation for the same sequence. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCliquePair.h>`_.

.. c:type:: AgnCliquePairClassification

  This enumerated type refers to all the possible outcomes when transcript cliques from two different sources are compared: ``AGN_CLIQUE_PAIR_UNCLASSIFIED``, ``AGN_CLIQUE_PAIR_PERFECT_MATCH``, ``AGN_CLIQUE_PAIR_MISLABELED``, ``AGN_CLIQUE_PAIR_CDS_MATCH``, ``AGN_CLIQUE_PAIR_EXON_MATCH``, ``AGN_CLIQUE_PAIR_UTR_MATCH``, and ``AGN_CLIQUE_PAIR_NON_MATCH``.



.. c:function:: void agn_clique_pair_build_model_vectors(AgnCliquePair *pair)

  Build a pair of model vectors to represent this pair of maximal transcripts or transcript cliques.

.. c:function:: AgnCliquePairClassification agn_clique_pair_classify(AgnCliquePair *pair)

  Based on the already-computed comparison statistics, classify this clique pair as a perfect match, a CDS match, etc. See :c:type:`AgnCliquePairClassification`.

.. c:function:: void agn_clique_pair_comparative_analysis(AgnCliquePair *pair)

  Compare the annotations for this pair, reference vs prediction.

.. c:function:: int agn_clique_pair_compare(void *p1, void *p2)

  Same as c:func:`agn_clique_pair_compare_direct`, but with pointer dereferencing.

.. c:function:: int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2)

  Determine which pair has higher comparison scores. Returns 1 if the first pair has better scores, -1 if the second pair has better scores, 0 if they are equal.

.. c:function:: int agn_clique_pair_compare_reverse(void *p1, void *p2)

  Negation of c:func:`agn_clique_pair_compare`.

.. c:function:: void agn_clique_pair_delete(AgnCliquePair *pair)

  Class destructor.

.. c:function:: double agn_clique_pair_get_edit_distance(AgnCliquePair *pair)

  Return the calculated annotation edit distance between this pair of transcripts or transcript cliques.

.. c:function:: AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair)

  Return a pointer to this pair's prediction transcript or transcript clique.

.. c:function:: const char *agn_clique_pair_get_pred_vector(AgnCliquePair *pair)

  Get the model vector associated with this pair's prediction transcript clique.

.. c:function:: AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair)

  Return a pointer to this pair's reference transcript or transcript clique.

.. c:function:: const char *agn_clique_pair_get_refr_vector(AgnCliquePair *pair)

  Get the model vector associated with this pair's reference transcript clique.

.. c:function:: AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair)

  Return a pointer to this clique pair's comparison statistics.

.. c:function:: bool agn_clique_pair_has_utrs(AgnCliquePair *pair)

  Determine whether there are UTRs in this clique pair

.. c:function:: bool agn_clique_pair_is_simple(AgnCliquePair *pair)

  Does this clique pair contain a single reference transcript and a single prediction transcript?

.. c:function:: GtUword agn_clique_pair_length(AgnCliquePair *pair)

  Get the length of the locus to which this clique pair belongs.

.. c:function:: bool agn_clique_pair_needs_comparison(AgnCliquePair *pair)

  Determine whether this clique pair needs comparison (i.e., whether there are both reference and prediction transcripts).

.. c:function:: AgnCliquePair* agn_clique_pair_new(const char *seqid, AgnTranscriptClique *refr_clique, AgnTranscriptClique *pred_clique, GtRange *locus_range)

  Class constructor.

.. c:function:: void agn_clique_pair_record_characteristics(AgnCliquePair *pair, AgnCompResultDesc *desc)

  Add information about a clique pair to a set of aggregate characteristics.

.. c:function:: bool agn_clique_pair_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Module AgnComparEval
--------------------

A collection of data structures and functions for comparative evaluation of annotations. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnComparEval.h>`_.

.. c:type:: AgnCompStatsScaled

  This struct is used to aggregate counts and statistics regarding the nucleotide-level comparison and analysis of gene structure.



.. c:type:: AgnCompStatsBinary

  This struct is used to aggregate counts and statistics regarding the structural-level comparison (i.e., at the level of whole CDS segments, whole exons, and whole UTRs) and analysis of gene structure.



.. c:type:: AgnCompSummary

  This struct contains various counts to be reported in the summary report.



.. c:type:: AgnComparison

  This struct aggregates all the counts and stats that go into a comparison, including structural-level and nucleotide-level counts and stats.



.. c:type:: AgnCompResultDesc

  Each transcript clique pair that is compared is classified as one of the following: perfect match; perfect match with mislabeled UTRs; CDS match; exon structure match; UTR structure match; non-match. When reporting the results of a comparative analysis, it may be useful to (as is done by ParsEval) show some basic information about clique pairs that fall under each classification category. The counts in this struct are necessary to calculate those summary characteristics.



.. c:type:: AgnCompResultSummary

  This struct is used to aggregate characteristics for all of the classification categories.



.. c:type:: AgnCompEvaluation

  This struct provides a convenient way to manage the counts, stats, and results corresponding to one or more comparisons.



.. c:type:: AgnCompareFilters

  This struct contains a list of filters to be used in determining which loci should be included/excluded in a comparative analysis.



.. c:function:: void agn_comp_evaluation_combine(AgnCompEvaluation *data, AgnCompEvaluation *data_to_add)

  Take values from one data set and add them to the other.

.. c:function:: void agn_comp_evaluation_init(AgnCompEvaluation *data)

  Initialize to default values.

.. c:function:: void agn_comp_result_summary_combine(AgnCompResultSummary *desc, AgnCompResultSummary *desc_to_add)

  Take values from one description and add them to the other.

.. c:function:: void agn_comp_result_summary_init(AgnCompResultSummary *desc)

  Initialize to default values.

.. c:function:: void agn_comp_result_desc_combine(AgnCompResultDesc *desc, AgnCompResultDesc *desc_to_add)

  Take the counts from one description and add them to a larger aggregate set of counts.

.. c:function:: void agn_comp_result_desc_init(AgnCompResultDesc *desc)

  Initialize characteristics to default values.

.. c:function:: void agn_comp_summary_combine(AgnCompSummary *s1, AgnCompSummary *s2)

  Take one set of values and add them to the other.

.. c:function:: void agn_comp_summary_init(AgnCompSummary *summary)

  Initialize to default values.

.. c:function:: void agn_comparison_combine(AgnComparison *c1, AgnComparison *c2)

  Take stats from one comparison and add them to the other.

.. c:function:: void agn_comparison_init(AgnComparison *comparison)

  Initialize comparison stats to default values.

.. c:function:: void agn_compare_filters_init(AgnCompareFilters *filters)

  Initialize filters to default values.

.. c:function:: void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream, AgnLogger *logger)

  Parse the filter configuration file (from ``instream``) to set the filters appropriately.

.. c:function:: void agn_comp_stats_binary_init(AgnCompStatsBinary *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats)

  Calculate stats from the given counts.

.. c:function:: void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats)

  Initialize comparison counts/stats to default values.

.. c:function:: void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats)

  Calculate stats from the given counts.

Class AgnFilterStream
---------------------

.. c:type:: AgnFilterStream

  Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream used to select features of a certain type from a node stream. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnFilterStream.h>`_.

.. c:function:: GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream, GtHashmap *typestokeep)

  Class constructor. The keys of the ``typestokeep`` hashmap should be the type(s) to be kept from the node stream. Any non-NULL value can be associated with those keys.

Class AgnGeneLocus
------------------

.. c:type:: AgnGeneLocus

  The purpose of the AgnGeneLocus class is to store all of the data associated with a distinct locus, in many cases to facilitate the comparison of two sets of gene structure annotations for that locus. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnGeneLocus.h>`_.

.. c:type:: AgnComparisonSource

  When tracking the source of an annotation for comparison purposes, use this enumerated type to refer to reference (``REFERENCESOURCE``) vs prediction (``PREDICTIONSOURCE``) annotations. ``DEFAULTSOURCE`` is for when the source is not a concern.



.. c:type:: AgnGeneLocusPngMetadata

  This data structure provides a convenient container for metadata needed to produce a PNG graphic for pairwise comparison loci.



.. c:type:: AgnGeneLocusSummary

  This data structure provides a summary of the data and comparisons associated with a given locus.



.. c:function:: void agn_gene_locus_add(AgnGeneLocus *locus, GtFeatureNode *gene, AgnComparisonSource source)

  Associate the given gene annotation with this gene locus. Rather than calling this function directly, users are recommended to use one of the following macros: ``agn_gene_locus_add_pred_gene(locus, gene)`` and ``agn_gene_locus_add_refr_gene(locus, gene)``, to be used when keeping track of an annotation's source is important (i.e. for pairwise comparison); and ``agn_gene_locus_add_gene(locus, gene)`` otherwise.

.. c:function:: void agn_gene_locus_aggregate_results(AgnGeneLocus *locus, AgnCompEvaluation *eval)

  Add the locus' comparison statistics to a set of aggregate statistics.

.. c:function:: AgnGeneLocus *agn_gene_locus_clone(AgnGeneLocus *locus)

  Do a semi-shallow copy of this data structure--for members whose data types support reference counting, the same pointer is used and the reference is incremented. For the other members a new object is created and populated with the same content.

.. c:function:: int agn_gene_locus_array_compare(const void *p1, const void *p2)

  Analog of ``strcmp`` for comparing AgnGeneLocus objects, used for sorting GtArray objects containing AgnGeneLocus objects.

.. c:function:: GtUword agn_gene_locus_cds_length(AgnGeneLocus *locus, AgnComparisonSource src)

  The combined length of all coding sequences associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_refr_cds_length(locus)`` for the combined length of all reference CDSs, ``agn_gene_locus_pred_cds_length(locus)`` for the combined length of all prediction CDSs, and ``agn_gene_locus_get_cds_length(locus)`` for the combined length of all CDSs.

.. c:function:: GtArray *agn_gene_locus_comparative_analysis(AgnGeneLocus *locus)

  Compare every reference transcript clique with every prediction transcript clique. For gene loci with multiple transcript cliques, each comparison is not necessarily reported. Instead, we report the set of clique pairs that provides the optimal pairing of reference and prediction transcripts. If there are more reference transcript cliques than prediction cliques (or vice versa), these unmatched cliques are reported separately.

.. c:function:: void agn_gene_locus_delete(AgnGeneLocus *locus)

  Class destructor.

.. c:function:: GtUword agn_gene_locus_exon_num(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the number of exons for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_num_pred_exons(locus)`` for the number of prediction exons, ``agn_gene_locus_num_refr_exons(locus)`` for the number of reference exons, or ``agn_gene_locus_num_exons(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: bool agn_gene_locus_filter(AgnGeneLocus *locus, AgnCompareFilters *filters)

  Given a set of filtering criteria, determine whether a locus meets those criteria. Returns true if the locus should be filtered (if it does not meet the criteria), false otherwise.

.. c:function:: GtArray *agn_gene_locus_genes(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the genes associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_pred_genes(locus)`` to retrieve prediction genes, ``agn_gene_locus_refr_genes(locus)`` to retrieve reference genes, or ``agn_gene_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src)

  Get IDs of the genes associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_pred_gene_ids(locus)`` to retrieve prediction genes IDs, ``agn_gene_locus_refr_gene_ids(locus)`` to retrieve reference genes IDs, or ``agn_gene_locus_get_gene_ids(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_gene_locus_gene_num(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the number of genes for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_num_pred_genes(locus)`` for the number of prediction genes, ``agn_gene_locus_num_refr_genes(locus)`` for the number of reference genes, or ``agn_gene_locus_num_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_gene_locus_get_end(AgnGeneLocus *locus)

  Get this locus' end coordinate.

.. c:function:: GtUword agn_gene_locus_get_length(AgnGeneLocus *locus)

  Get this locus' length.

.. c:function:: const char* agn_gene_locus_get_seqid(AgnGeneLocus *locus)

  Get this locus' sequence ID.

.. c:function:: GtUword agn_gene_locus_get_start(AgnGeneLocus *locus)

  Get this locus' start coordinate.

.. c:function:: GtArray *agn_gene_locus_get_unique_pred_cliques(AgnGeneLocus *locus)

  Get a list of all the prediction transcript cliques that have no corresponding reference transcript clique.

.. c:function:: GtArray *agn_gene_locus_get_unique_refr_cliques(AgnGeneLocus *locus)

  Get a list of all the reference transcript cliques that have no corresponding prediction transcript clique.

.. c:function:: AgnGeneLocus* agn_gene_locus_new(const char *seqid)

  Class constructor.

.. c:function:: GtUword agn_gene_locus_num_clique_pairs(AgnGeneLocus *locus)

  Report the number of clique pairs to be reported for this locus.

.. c:function:: void agn_gene_locus_png_track_selector(GtBlock *block, GtStr *track,void *data)

  Track selector function for generating PNG graphics of pairwise comparison loci. The track name to will be written to ``track``.

.. c:function:: void agn_gene_locus_print_png(AgnGeneLocus *locus, AgnGeneLocusPngMetadata *metadata)

  Print a PNG graphic for this locus.

.. c:function:: GtRange agn_gene_locus_range(AgnGeneLocus *locus)

  Return the coordinates of this locus.

.. c:function:: void agn_gene_locus_set_range(AgnGeneLocus *locus, GtUword start, GtUword end)

  Set the range of this locus, no questions asked.

.. c:function:: double agn_gene_locus_splice_complexity(AgnGeneLocus *locus, AgnComparisonSource src)

  Calculate the splice complexity of this gene locus. Rather than calling this method directly, users are recommended to use one of the following macros: ``agn_gene_locus_prep_splice_complexity(locus)`` to calculate the splice complexity of just the prediction transcripts, ``agn_gene_locus_refr_splice_complexity(locus)`` to calculate the splice complexity of just the reference transcripts, and ``agn_gene_locus_calc_splice_complexity(locus)`` to calculate the splice complexity taking into account all transcripts.

.. c:function:: void agn_gene_locus_summary_init(AgnGeneLocusSummary *summary)

  Class constructor.

.. c:function:: void agn_gene_locus_to_gff3(AgnGeneLocus *locus, FILE *outstream, const char *source)

  Print the locus in GFF3 format. If ``source`` is NULL, the string "AEGeAn" will be used.

.. c:function:: GtArray *agn_gene_locus_transcripts(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the transcripts associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_pred_transcripts(locus)`` to retrieve prediction transcripts, ``agn_gene_locus_refr_transcripts(locus)`` to retrieve reference transcripts, or ``agn_gene_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_gene_locus_transcript_ids(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the transcript IDs associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_gene_locus_pred_transcripts(locus)`` to retrieve prediction IDs, ``agn_gene_locus_refr_transcripts(locus)`` to retrieve reference IDs, or ``agn_gene_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_gene_locus_transcript_num(AgnGeneLocus *locus, AgnComparisonSource src)

  Get the number of transcripts for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_transcript_locus_num_pred_transcripts(locus)`` for the number of prediction transcripts, ``agn_transcript_locus_num_refr_transcripts(locus)`` for the number of reference transcripts, or ``agn_transcript_locus_num_transcripts(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: bool agn_gene_locus_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Module AgnGtExtensions
----------------------

A collection of extensions to core GenomeTools classes. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnGtExtensions.h>`_.

.. c:function:: GtArray* agn_gt_array_copy(GtArray *source, size_t size)

  This function makes a copy of an array.

.. c:function:: void agn_gt_feature_index_to_gff3(GtFeatureIndex *index, FILE *outstream)

  Write the given feature index to GFF3 format.

.. c:function:: GtUword agn_gt_feature_node_cds_length(GtFeatureNode *transcript)

  Calculate the length of the given transcript's coding sequence in amino acids.

.. c:function:: GtArray *agn_gt_feature_node_children_of_type(GtFeatureNode *fn, bool (*typetestfunc)(GtFeatureNode *))

  Gather the children of a given feature that have a certain type. Type is tested by ``typetestfunc``, which accepts a single ``GtFeatureNode`` object.

.. c:function:: bool agn_gt_feature_node_fix_parent_attribute(GtFeatureNode *feature, GtFeatureNode *parent)

  When a feature has multiple parents but only one of them is valid, this function will fix the ``Parent`` attribute so that it only points to the valid parent.

.. c:function:: bool agn_gt_feature_node_is_cds_feature(GtFeatureNode *fn)

  Determine whether the given feature belongs to a CDS.

.. c:function:: bool agn_gt_feature_node_is_exon_feature(GtFeatureNode *fn)

  Determine whether the given feature is an exon.

.. c:function:: bool agn_gt_feature_node_is_gene_feature(GtFeatureNode *fn)

  Determine whether the given feature is a gene.

.. c:function:: bool agn_gt_feature_node_is_intron_feature(GtFeatureNode *fn)

  Determine whether the given feature is an intron.

.. c:function:: bool agn_gt_feature_node_is_mrna_feature(GtFeatureNode *fn)

  Determine whether the given feature is an mRNA.

.. c:function:: bool agn_gt_feature_node_is_start_codon_feature(GtFeatureNode *fn)

  Determine whether the given feature is a start codon.

.. c:function:: bool agn_gt_feature_node_is_stop_codon_feature(GtFeatureNode *fn)

  Determine whether the given feature is a stop codon.

.. c:function:: bool agn_gt_feature_node_is_utr_feature(GtFeatureNode *fn)

  Determine whether the given feature is part of a UTR.

.. c:function:: GtUword agn_gt_feature_node_num_transcripts(GtFeatureNode *gene)

  Determine the number of transcripts for the given gene feature.

.. c:function:: bool agn_gt_feature_node_overlap(GtFeatureNode *first, GtFeatureNode *second)

  Determine whether the given features overlap.

.. c:function:: bool agn_gt_feature_node_range_contains(GtFeatureNode *n1, GtFeatureNode *n2)

  Determine whether the range of n2 falls within the range of n1.

.. c:function:: void agn_gt_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)

  The ``gt_feature_node_remove_leaf`` function only allows removal of a leaf node. This function will remove all of a node's children so that it is a leaf node, which can then be removed.

.. c:function:: void agn_gt_feature_node_set_source_recursive(GtFeatureNode *feature, GtStr *source)

  Reset the source of the feature and all its children to the given value.

.. c:function:: void agn_gt_feature_node_to_gff3(GtFeatureNode *feature, FILE *outstream, bool printchildren, char *prefix, GtHashmap *filtered_types)

  Print the given feature in GFF3 format. Use ``filtered_types`` to provide a list of feature types to exclude when printing. If ``filtered_types`` is NULL, a default exclusion list will be used.

.. c:function:: int agn_gt_genome_node_compare(const void *n1, const void *n2)

  Comparison function to be used for sorting GtGenomeNode objects stored in a GtArray (for GtDlist, use gt_genome_node_cmp).

.. c:function:: char agn_gt_phase_to_char(GtPhase phase)

  Convert a GtPhase object into its corresponding character representation.

.. c:function:: char agn_gt_strand_to_char(GtStrand strand)

  Convert a GtStrand object into its corresponding character representation.

.. c:function:: GtStrArray* agn_gt_str_array_intersection(GtStrArray *a1, GtStrArray *a2)

  Find the strings that are common to both string arrays.

.. c:function:: GtStrArray* agn_gt_str_array_union(GtStrArray *a1, GtStrArray *a2)

  Find the strings that are present in either (or both) of the string arrays.

Class AgnInferCDSVisitor
------------------------

.. c:type:: AgnInferCDSVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring an mRNA's CDS from explicitly defined exon and start/stop codon features. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferCDSVisitor.h>`_.

.. c:function:: GtNodeVisitor* agn_infer_cds_visitor_new(AgnLogger *logger)

  Class constructor.

.. c:function:: bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class.

Class AgnInferExonsVisitor
--------------------------

.. c:type:: AgnInferExonsVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring exon features when only CDS and UTR features are provided explicitly.  See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferExonsVisitor.h>`_.

.. c:function:: GtNodeVisitor* agn_infer_exons_visitor_new(AgnLogger *logger)

  Class constructor.

.. c:function:: bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class.

Class AgnLocusIndex
-------------------

.. c:type:: AgnLocusIndex

  FIXME See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocusIndex.h>`_.

.. c:type:: typedef void (*AgnLocusIndexVisitFunc)(AgnGeneLocus *, void *)

  Signature functions must match to be applied to each locus in the index. The function will be called once for each locus, which will be passed as the first argument to the function. a second argument is available for an optional pointer to supplementary data (if needed). See :c:func:`agn_locus_index_comparative_analysis`.

.. c:function:: void agn_locus_index_comparative_analysis(AgnLocusIndex *idx, const char *seqid, AgnLocusIndexVisitFunc preanalyfunc, AgnLocusIndexVisitFunc postanalyfunc, void *analyfuncdata, AgnLogger *logger)

  Perform a comparative analysis of each locus associated with ``seqid`` in this index. If ``preanalyfunc`` is not NULL, it will be applied to each locus immediately before comparative analysis. If ``postanalyfunc`` is not NULL, it will be applied to each locus immediately following comparative analysis. ``analyfuncdata`` will be passed as supplementary data to both functions.

.. c:function:: void agn_locus_index_delete(AgnLocusIndex *idx)

  Class destructor.

.. c:type:: void agn_locus_index_find(AgnLocusIndex *idx, const char *seqid, GtRange *range, GtArray *loci)

  Find all overlapping features in the given range stored in this locus index and store them in ``loci``.

.. c:function:: GtArray *agn_locus_index_get(AgnLocusIndex *idx, const char *seqid)

  Retrieve all loci corresponding to the specified sequence ID.

.. c:type:: GtArray *agn_locus_index_interval_loci(AgnLocusIndex *idx, const char *seqid, GtUword delta, bool skipterminal)

  Compute interval loci with the given ``delta``. If running on incomplete (contig/scaffold) genomic sequences, consider setting ``skipterminal`` to true to ignore the ends of the sequence.

.. c:function:: AgnLocusIndex *agn_locus_index_new(bool freeondelete)

  Class constructor

.. c:function:: GtUword agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx, GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, AgnCompareFilters *filters, AgnLogger *logger)

  Given a pair of annotation feature sets in memory, identify loci while keeping the two sources of annotation separate (to enable comparison).

.. c:function:: GtUword agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx, const char *refrfile, const char *predfile, AgnCompareFilters *filters, AgnLogger *logger)

  Given a pair of annotation feature sets in memory, identify loci while keeping the two sources of annotation separate (to enable comparison).

.. c:function:: GtUword agn_locus_index_parse_memory(AgnLocusIndex *idx, GtFeatureIndex *features, AgnLogger *logger)

  Identify loci given an index of annotation features.

.. c:function:: GtUword agn_locus_index_parse_disk(AgnLocusIndex *idx, int numfiles, const char **filenames, AgnLogger *logger)

  Identify loci from the given set of annotation files.

.. c:function:: GtStrArray *agn_locus_index_seqids(AgnLocusIndex *idx)

  Get a list of the seqids stored in this locus index.

Class AgnLogger
---------------

.. c:type:: AgnLogger

  The AgnLogger class is desiged to store error, warning, and status messages. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLogger.h>`_.

.. c:function:: void agn_logger_delete(AgnLogger *logger)

  Class destructor.

.. c:function:: GtArray *agn_logger_get_error_messages(AgnLogger *logger)

  Return an array containing all the error messages associated with this logger. The user is responsible for freeing this array, but not its contents.

.. c:function:: GtArray *agn_logger_get_status_messages(AgnLogger *logger)

  Return an array containing all the status messages associated with this logger. The user is responsible for freeing this array, but not its contents.

.. c:function:: GtArray *agn_logger_get_warning_messages(AgnLogger *logger)

  Return an array containing all the warning messages associated with this logger. The user is responsible for freeing this array, but not its contents.

.. c:function:: bool agn_logger_has_error(AgnLogger *logger)

  Have any errors been logged?

.. c:function:: bool agn_logger_has_status(AgnLogger *logger)

  Have any status messages been logged?

.. c:function:: bool agn_logger_has_warning(AgnLogger *logger)

  Have any warnings been logged?

.. c:function:: void agn_logger_log_error(AgnLogger *logger, const char *format, ...)

  Add an error message to the logger using ``printf``-style string formatting.

.. c:function:: void agn_logger_log_status(AgnLogger *logger, const char *format, ...)

  Add a status message/update to the logger using ``printf``-style string formatting.

.. c:function:: void agn_logger_log_warning(AgnLogger *logger, const char *format, ...)

  Add a warning message to the logger using ``print``-style string formatting.

.. c:function:: AgnLogger *agn_logger_new()

  Class constructor.

.. c:function:: bool agn_logger_print_all(AgnLogger *logger, FILE *outstream, const char *format, ...)

  Print the status messages, warnings, and errors that have been logged to the given file stream, ``printf``-style. Returns true if any errors were printed, false otherwise.

.. c:function:: bool agn_logger_print_error(AgnLogger *logger, FILE *outstream, const char *format, ...)

  Print the error messages associated with this logger to the given file stream. Returns true if errors were printed.

.. c:function:: bool agn_logger_print_status(AgnLogger *logger, FILE *outstream, const char *format, ...)

  Print the status messages associated with this logger to the given file stream. Returns true if any messages were printed, false otherwise.

.. c:function:: bool agn_logger_print_warning(AgnLogger *logger, FILE *outstream, const char *format, ...)

  Print the warning messages associated with this logger to the given file stream. Returns true if any warnings were printed, false otherwise.

.. c:function:: void agn_logger_unset(AgnLogger *logger)

  Reset this logger object.

Module AgnTestData
------------------

A collection of functions facilitating unit testing of various AEGeAn classes and modules. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTestData.h>`_.

.. c:function:: GtArray *agn_test_data_grape()

  Example from grape.

.. c:function:: GtArray *agn_test_data_grape_codons()

  Example from grape: gene structure annotated with exon and start / stop codon features--CDS is implicitly defined by these features.

.. c:function:: GtArray *agn_test_data_grape_sansexons()

  Example from grape: gene structure annotated with CDS and UTR features--exons are implicitly defined by these features.

.. c:function:: GtFeatureNode *agn_test_data_eden()

  Create the canonical gene structure (from the GFF3 specification) in memory.

Class AgnTranscriptClique
-------------------------

.. c:type:: AgnTranscriptClique

  The purpose of the AgnTranscriptClique class is to store data pertaining to an individual maximal transcript clique. This clique may only contain a single transcript, or it may contain many. The only stipulation is that the transcripts do not overlap. See the `class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTranscriptClique.h>`_.

.. c:type:: typedef void (*AgnCliqueVisitFunc)(GtFeatureNode*, void*)

   The signature that functions must match to be applied to each transcript in the given clique. The function will be called once for each transcript in the clique. The transcript will be passed as the first argument, and a second argument is available for an optional pointer to supplementary data (if needed). See :c:func:`agn_transcript_clique_traverse`.

.. c:function:: void agn_transcript_clique_add(AgnTranscriptClique *clique, GtFeatureNode *transcript)

  Add a transcript to this clique.

.. c:function:: GtUword agn_transcript_clique_cds_length(AgnTranscriptClique *clique)

  Get the CDS length (in amino acids) for this transcript clique.

.. c:function:: AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)

  Make a shallow copy of this transcript clique.

.. c:function:: void agn_transcript_clique_delete(AgnTranscriptClique *clique)

  Class destructor.

.. c:function:: bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique, GtHashmap *map)

  Determine whether any of the transcript IDs associated with this clique are keys in the given hash map.

.. c:function:: const char *agn_transcript_clique_id(AgnTranscriptClique *clique)

  Retrieve the ID attribute of the transcript associated with this clique. Will cause an assertion error if there is more than one trancript associated with the clique.

.. c:function:: AgnTranscriptClique* agn_transcript_clique_new()

  Class constructor.

.. c:function:: GtUword agn_transcript_clique_num_exons(AgnTranscriptClique *clique)

  Get the number of exons in this clique.

.. c:function:: GtUword agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)

  Get the number of UTR segments in this clique.

.. c:function:: void agn_transcript_clique_print_ids(AgnTranscriptClique *clique, FILE *outstream)

  Print the IDs of this clique's transcripts to the given output stream.

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

Module AgnUtils
---------------

A collection of assorted core utility functions. See the `module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUtils.h>`_.

.. c:type:: AgnSequenceRegion

  Simple data structure for referencing genomic locations.



.. c:function:: void agn_bron_kerbosch( GtArray *R, GtArray *P, GtArray *X, GtArray *cliques, bool skipsimplecliques )

  The Bron-Kerbosch algorithm is an algorithm for enumerating all maximal cliques in an undirected graph. See the `algorithm's Wikipedia entry <http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm>`_ for a description of ``R``, ``P``, and ``X``. All maximal cliques will be stored in ``cliques``. If ``skipsimplecliques`` is true, cliques containing a single item will not be stored.

.. c:function:: double agn_calc_edit_distance(GtFeatureNode *t1, GtFeatureNode *t2)

  A paper by Eilbeck `et al` (http://dx.doi.org/10.1186/1471-2105-10-67) described the annotation edit distance as a measure for comparative evaluation of annotations. This function calculates the AED between the two given annotations.

.. c:function:: double agn_calc_splice_complexity(GtArray *transcripts)

  Determine the splice complexity of the given set of transcripts.

.. c:function:: GtArray* agn_enumerate_feature_cliques(GtArray *feature_set)

  If reference transcripts belonging to the same locus overlap, they must be separated before comparison with prediction transcript models (and vice versa). This is an instance of the maximal clique enumeration problem (NP-complete), for which the Bron-Kerbosch algorithm provides a solution.

.. c:function:: GtArray* agn_feature_neighbors(GtGenomeNode *feature, GtArray *feature_set)

  For a set of features, we can construct a graph where each node represents a feature and where two nodes are connected if the corresponding features do not overlap. This function returns the intersection of feature_set with the neighbors of feature (where a "neighbor" refers to an adjacent node).

.. c:function:: FILE *agn_fopen(const char *filename, const char *mode, FILE *errstream)

  Wrapper around the stdio.h function that will exit in case of an IO error.

.. c:function:: GtFeatureIndex *agn_import_canonical(int numfiles, const char **filenames, AgnLogger *logger)

  Load canonical protein-coding genes from the given GFF3 files into memory.

.. c:function:: GtFeatureIndex *agn_import_simple(int numfiles, const char **filenames, char *type, AgnLogger *logger)

  Load features whose type is equal to ``type`` into memory from the given GFF3 files.

.. c:function:: bool agn_infer_cds_range_from_exon_and_codons(GtRange *exon_range, GtRange *leftcodon_range, GtRange *rightcodon_range, GtRange *cds_range)

  Given an exon and the start/stop codons associated with its corresponding mRNA, determine which parts of the exon (if any) correspond to coding sequence. If the exon contains coding sequence, the range of that coding sequence will be stored in ``cds_range`` and the function will return true. Otherwise, the function will return false. If the mRNA is on the forward strand, ``left_codon_range`` should contain the coordinates for the start codon and ``right_codon_range`` should contain coordinates for the stop codon. If the mRNA is on the reverse strand, these should be swapped.

.. c:function:: GtStrArray* agn_seq_intersection(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, AgnLogger *logger)

  Given two feature indices, determine which sequences are common between them and return those sequences' IDs as a string array.

.. c:function:: GtStrArray* agn_seq_union(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, AgnLogger *logger)

  Given two feature indices, determine all of the sequences that are annotated by either of them and return those sequences' IDs as a string array.

.. c:function:: int agn_sprintf_comma(GtUword n, char *buffer)

  Format the given non-negative number with commas as the thousands separator. The resulting string will be written to ``buffer``.

.. c:function:: int agn_string_compare(const void *p1, const void *p2)

  Dereference the given pointers and compare the resulting strings (a la ``strcmp``).

.. c:function:: GtRange agn_transcript_cds_range(GtFeatureNode *transcript)

  Determine the start and end coordinates of the given transcript's CDS.

.. c:function:: void agn_transcript_structure_gbk(GtFeatureNode *transcript, FILE *outstream)

  Write the structure of a gene transcript in GenBank format to ``outstream``.

