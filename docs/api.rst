AEGeAn C API
============

The AEGeAn Toolkit relies heavily on data
types implemented by the GenomeTools library. For data types beginning with
``Gt``, see the GenomeTools API documentation at
http://genometools.org/libgenometools.html.

Class AgnCliquePair
-------------------

.. c:type:: AgnCliquePair

  The AgnCliquePair class facilitates comparison of two alternative sources of annotation for the same sequence. See the `AgnCliquePair class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCliquePair.h>`_.

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

.. c:function:: AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair)

  Return a pointer to this clique pairs comparison statistics.

.. c:function:: AgnCliquePair* agn_clique_pair_new(AgnTranscriptClique *refr, AgnTranscriptClique *pred)

  Class constructor.

.. c:function:: bool agn_clique_pair_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnCompareReportHTML
--------------------------

.. c:type:: AgnCompareReportHTML

  The ``AgnCompareReportHTML`` class is an extension of the ``AgnCompareReport`` class. This node visitor relies on its parent class to process a stream of ``AgnLocus`` objects (containing two alternative sources of annotation to be compared) and then produces textual reports of the comparison statistics. See the `AgnCompareReportHTML class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCompareReportHTML.h>`_.

.. c:type:: typedef void (*AgnCompareReportHTMLOverviewFunc)(FILE *outstream, void *data)

  By default, the ParsEval summary report includes an overview with the start time, filenames, and command-line arguments. Users can override this behavior by specifying a callback function that follows this signature.

.. c:function:: void agn_compare_report_html_create_summary(AgnCompareReportHTML *rpt)

  After the node stream has been processed, call this function to write a summary of all locus comparisons to the output directory.

.. c:function:: GtNodeVisitor *agn_compare_report_html_new(const char *outdir, bool gff3, AgnLocusPngMetadata *pngdata, GtLogger *logger)

  Class constructor. Creates a node visitor used to process a stream of ``AgnLocus`` objects containing two sources of annotation to be compared. Reports will be written in ``outdir`` and status messages will be written to the logger.

.. c:function:: void agn_compare_report_html_reset_summary_title(AgnCompareReportHTML *rpt, GtStr *title_string)

  By default, the summary report's title will be 'ParsEval Summary'. Use this function to replace the title text.

.. c:function:: void agn_compare_report_html_set_overview_func(AgnCompareReportHTML *rpt, AgnCompareReportHTMLOverviewFunc func, void *funcdata)

  Specify a callback function to be used when printing an overview on the summary report.

Class AgnCompareReportText
--------------------------

.. c:type:: AgnCompareReportText

  The ``AgnCompareReportText`` class is an extension of the ``AgnCompareReport`` class. This node visitor relies on its parent class to process a stream of ``AgnLocus`` objects (containing two alternative sources of annotation to be compared) and then produces textual reports of the comparison statistics. See the `AgnCompareReportText class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCompareReportText.h>`_.

.. c:function:: void agn_compare_report_text_create_summary(AgnCompareReportText *rpt, FILE *outstream)

  After the node stream has been processed, call this function to write a summary of all locus comparisons to ``outstream``.

.. c:function:: GtNodeVisitor *agn_compare_report_text_new(FILE *outstream, bool gff3, GtLogger *logger)

  Class constructor. Creates a node visitor used to process a stream of ``AgnLocus`` objects containing two sources of annotation to be compared. Reports will be written to ``outstream`` and status messages will be written to the logger.

Module AgnComparison
--------------------

Data structures and functions related to comparative assessment of gene/transcript annotations. See the `AgnComparison module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnComparison.h>`_.

.. c:type:: AgnCompStatsBinary

  This struct is used to aggregate counts and statistics regarding the structural-level comparison (i.e., at the level of whole CDS segments, whole exons, and whole UTRs) and analysis of gene structure. See header file for details.



.. c:type:: AgnCompStatsScaled

  This struct is used to aggregate counts and statistics regarding the nucleotide-level comparison and analysis of gene structure. See header file for details.



.. c:type:: AgnComparison

  This struct aggregates all the counts and stats that go into a comparison, including structural-level and nucleotide-level counts and stats. See header file for details.



.. c:type:: AgnCompClassification

  This enumerated type refers to all the possible outcomes when annotations from two different sources are compared: ``AGN_COMP_CLASS_UNCLASSIFIED``, ``AGN_COMP_CLASS_PERFECT_MATCH``, ``AGN_COMP_CLASS_MISLABELED``, ``AGN_COMP_CLASS_CDS_MATCH``, ``AGN_COMP_CLASS_EXON_MATCH``, ``AGN_COMP_CLASS_UTR_MATCH``, and ``AGN_COMP_CLASS_NON_MATCH``.



.. c:type:: AgnCompInfo

  This struct contains various counts to be reported in the summary report.



.. c:type:: AgnCompClassDesc

  When reporting the results of a comparative analysis, it may be useful to (as is done by ParsEval) show some basic information about clique pairs that fall under each classification category. The counts in this struct are necessary to calculate those summary characteristics.



.. c:type:: AgnCompClassSummary

  This struct is used to aggregate descriptions for all of the classification categories.



.. c:type:: AgnComparisonData

  Aggregate various data related to comparison of annotations.



.. c:function:: void agn_comparison_aggregate(AgnComparison *agg_cmp, AgnComparison *cmp)

  Function used to combine similarity stats from many different comparisons into a single aggregate summary.

.. c:function:: void agn_comparison_data_aggregate(AgnComparisonData *agg_data, AgnComparisonData *data)

  Add counts and stats from ``data`` to ``agg_data``.

.. c:function:: void agn_comparison_data_init(AgnComparisonData *data)

  Initialize counts and stats to default values.

.. c:function:: void agn_comparison_init(AgnComparison *comparison)

  Initialize comparison stats to default values.

.. c:function:: void agn_comparison_print(AgnComparison *stats, FILE *outstream)

  Print the comparison stats to the given file.

.. c:function:: void agn_comparison_resolve(AgnComparison *comparison)

  Calculate stats from the given counts.

.. c:function:: bool agn_comparison_test(AgnComparison *c1, AgnComparison *c2)

  Returns true if c1 and c2 contain identical values, false otherwise.

.. c:function:: void agn_comp_class_desc_aggregate(AgnCompClassDesc *agg_desc, AgnCompClassDesc *desc)

  Add values from ``desc`` to ``agg_desc``.

.. c:function:: void agn_comp_class_desc_init(AgnCompClassDesc *desc)

  Initialize to default values.

.. c:function:: void agn_comp_class_summary_aggregate(AgnCompClassSummary *agg_summ, AgnCompClassSummary *summ)

  Add values from ``summ`` to ``agg_summ``.

.. c:function:: void agn_comp_class_summary_init(AgnCompClassSummary *summ)

  Initialize to default values.

.. c:function:: void agn_comp_info_aggregate(AgnCompInfo *agg_info, AgnCompInfo *info)

  Add values from ``info`` to ``agg_info``.

.. c:function:: void agn_comp_info_init(AgnCompInfo *info)

  Initialize to default values.

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

  Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream used to select features of a certain type from a node stream. See the `AgnFilterStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnFilterStream.h>`_.

.. c:function:: GtNodeStream* agn_filter_stream_new(GtNodeStream *in_stream, GtHashmap *typestokeep)

  Class constructor. The keys of the ``typestokeep`` hashmap should be the type(s) to be kept from the node stream. Any non-NULL value can be associated with those keys.

.. c:function:: bool agn_filter_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnGeneStream
-------------------

.. c:type:: AgnGeneStream

  Implements the ``GtNodeStream`` interface. Searches the complete feature graph of each feature node in the input for canonical protein-coding gene features. Some basic sanity checks are performed on the mRNA(s) associated with each gene, and genes are only delivered to the output stream if they include one or more valid mRNA subfeatures. See the `AgnGeneStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnGeneStream.h>`_.

.. c:function:: GtNodeStream* agn_gene_stream_new(GtNodeStream *in_stream, GtLogger *logger)

  Class constructor.

.. c:function:: bool agn_gene_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnInferCDSVisitor
------------------------

.. c:type:: AgnInferCDSVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring an mRNA's CDS from explicitly defined exon and start/stop codon features. See the `AgnInferCDSVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferCDSVisitor.h>`_.

.. c:function:: GtNodeStream* agn_infer_cds_stream_new(GtNodeStream *in, GtLogger *logger)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_infer_cds_visitor_new(GtLogger *logger)

  Constructor for the node visitor.

.. c:function:: bool agn_infer_cds_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnInferExonsVisitor
--------------------------

.. c:type:: AgnInferExonsVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for inferring exon features when only CDS and UTR features are provided explicitly.  See the `AgnInferExonsVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferExonsVisitor.h>`_.

.. c:function:: GtNodeStream* agn_infer_exons_stream_new(GtNodeStream *in, GtLogger *logger)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor* agn_infer_exons_visitor_new(GtLogger *logger)

  Class constructor for the node visitor.

.. c:function:: bool agn_infer_exons_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class.

Class AgnInferParentStream
--------------------------

.. c:type:: AgnInferParentStream

  Implements the GenomeTools ``GtNodeStream`` interface. This node stream creates new features as parents for the specified types. For example, if ``type_parents`` includes an entry with ``tRNA`` as the key and ``gene`` as the value, this node stream will create a ``gene`` feature for any ``tRNA`` feature that lacks a gene parent. See the `AgnInferParentStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferParentStream.h>`_.

.. c:function:: GtNodeStream* agn_infer_parent_stream_new(GtNodeStream *in_stream, GtHashmap *type_parents)

  Class constructor. The hashmap contains a list of key-value pairs, both strings. Any time the stream encounters a top-level (parentless) feature whose type is a key in the hashmap, a parent will be created for this feature of the type associated with the key.

.. c:function:: bool agn_infer_parent_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnIntervalLocusStream
----------------------------

.. c:type:: AgnIntervalLocusStream

  Implements the ``GtNodeStream`` interface. Input is a stream of gene/transcript loci and output is a stream of interval loci. See online docs for more information about interval loci (iLoci). See the `AgnIntervalLocusStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnIntervalLocusStream.h>`_.

.. c:function:: GtNodeStream *agn_interval_locus_stream_new(GtNodeStream *locus_stream, GtUword delta, int endmode, GtLogger *logger)

  Class constructor. The delta parameter specifies how far beyond each transcript the iLocus boundaries should extend, and the minimum length of an iLocus containing no transcripts. If ``endmode == 0``, all iLoci will be included in the output; if ``endmode < 0``, terminal iLoci will not be included in the output; and if ``endmode > 0``, then only terminal iLoci will be included in the output. See the online docs for a complete description of iLoci.

.. c:function:: void agn_interval_locus_stream_set_idformat(AgnIntervalLocusStream *stream, const char *format)

  iLoci created by this stream are assigned an ID with an arbitrary number. The default format is 'iLocus%lu' (that is, iLocus1, iLocus2, etc). Use this function to override the default ID format.

.. c:function:: void agn_interval_locus_stream_set_source(AgnIntervalLocusStream *stream, GtStr *source)

  Set the source value to be used for all iLoci created by this stream. Default value is 'AEGeAn::AgnIntervalLocusStream'.

.. c:function:: bool agn_interval_locus_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnLocus
--------------

.. c:type:: AgnLocus

  The AgnLocus class represents gene loci and interval loci in memory and can be used to facilitate comparison of two different sources of annotation. Under the hood, each ``AgnLocus`` object is a feature node with one or more gene features as direct children. See the `AgnLocus class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocus.h>`_.

.. c:type:: AgnComparisonSource

  When tracking the source of an annotation for comparison purposes, use this enumerated type to refer to reference (``REFERENCESOURCE``) vs prediction (``PREDICTIONSOURCE``) annotations. ``DEFAULTSOURCE`` is for when the source is not a concern.



.. c:type:: AgnLocusPngMetadata

  This data structure provides a convenient container for metadata needed to produce a PNG graphic for pairwise comparison loci.



.. c:type:: AgnLocusFilterOp

  Comparison operators to use when filtering loci.



.. c:type:: AgnLocusFilter

  Data by which to filter a locus. If the value returned by ``function`` satisfies the criterion specified by ``testvalue`` and ``operator``, then the locus is to be kept.



.. c:function:: void agn_locus_add(AgnLocus *locus, GtFeatureNode *feature, AgnComparisonSource source)

  Associate the given annotation with this locus. Rather than calling this function directly, users are recommended to use one of the following macros: ``agn_locus_add_pred_feature(locus, gene)`` and ``agn_locus_add_refr_feature(locus, gene)``, to be used when keeping track of an annotation's source is important (i.e. for pairwise comparison); and ``agn_locus_add_feature(locus, gene)`` otherwise.

.. c:function:: AgnLocus *agn_locus_clone(AgnLocus *locus)

  Do a semi-shallow copy of this data structure--for members whose data types support reference counting, the same pointer is used and the reference is incremented. For the other members a new object is created and populated with the same content.

.. c:function:: GtUword agn_locus_cds_length(AgnLocus *locus, AgnComparisonSource src)

  The combined length of all coding sequences associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_refr_cds_length(locus)`` for the combined length of all reference CDSs, ``agn_locus_pred_cds_length(locus)`` for the combined length of all prediction CDSs, and ``agn_locus_get_cds_length(locus)`` for the combined length of all CDSs.

.. c:function:: void agn_locus_comparative_analysis(AgnLocus *locus, GtLogger *logger)

  Compare every reference transcript clique with every prediction transcript clique. For gene loci with multiple transcript cliques, each comparison is not necessarily reported. Instead, we report the set of clique pairs that provides the optimal pairing of reference and prediction transcripts. If there are more reference transcript cliques than prediction cliques (or vice versa), these unmatched cliques are reported separately.

.. c:function:: int agn_locus_array_compare(const void *p1, const void *p2)

  Analog of ``strcmp`` for sorting AgnLocus objects. Loci are first sorted lexicographically by sequence ID, and then spatially by genomic coordinates.

.. c:function:: void agn_locus_comparison_aggregate(AgnLocus *locus, AgnComparison *comp)

  Add this locus' internal comparison stats to a larger set of aggregate stats.

.. c:function:: void agn_locus_data_aggregate(AgnLocus *locus, AgnComparisonData *data)

  Add this locus' internal comparison stats to a larger set of aggregate stats.

.. c:function:: void agn_locus_delete(AgnLocus *locus)

  Class destructor.

.. c:function:: GtUword agn_locus_exon_num(AgnLocus *locus, AgnComparisonSource src)

  Get the number of exons for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_num_pred_exons(locus)`` for the number of prediction exons, ``agn_locus_num_refr_exons(locus)`` for the number of reference exons, or ``agn_locus_num_exons(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: void agn_locus_filter_parse(FILE *filterfile, GtArray *filters)

  Parse filters from ``filterfile`` and place ``AgnLocusFilter`` objects in ``filters``.

.. c:function:: bool agn_locus_filter_test(AgnLocus *locus, AgnLocusFilter *filter)

  Return true if ``locus`` satisfies the given filtering criterion.

.. c:function:: GtArray *agn_locus_get_unique_pred_cliques(AgnLocus *locus)

  Get a list of all the prediction transcript cliques that have no corresponding reference transcript clique.

.. c:function:: GtArray *agn_locus_get_unique_refr_cliques(AgnLocus *locus)

  Get a list of all the reference transcript cliques that have no corresponding prediction transcript clique.

.. c:function:: GtArray *agn_locus_genes(AgnLocus *locus, AgnComparisonSource src)

  Get the genes associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_genes(locus)`` to retrieve prediction genes, ``agn_locus_refr_genes(locus)`` to retrieve reference genes, or ``agn_locus_get_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_locus_gene_ids(AgnLocus *locus, AgnComparisonSource src)

  Get the gene IDs associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_gene_ids(locus)`` to retrieve prediction IDs, ``agn_locus_refr_gene_ids(locus)`` to retrieve reference IDs, or ``agn_locus_get_gene_ids(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_locus_gene_num(AgnLocus *locus, AgnComparisonSource src)

  Get the number of genes for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_num_pred_genes(locus)`` for the number of prediction genes, ``agn_locus_num_refr_genes(locus)`` for the number of reference genes, or ``agn_locus_num_genes(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_locus_mrnas(AgnLocus *locus, AgnComparisonSource src)

  Get the mRNAs associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_mrnas(locus)`` to retrieve prediction mRNAs, ``agn_locus_refr_mrnas(locus)`` to retrieve reference mRNAs, or ``agn_locus_get_mrnas(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtArray *agn_locus_mrna_ids(AgnLocus *locus, AgnComparisonSource src)

  Get the mRNA IDs associated with this locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_pred_mrna_ids(locus)`` to retrieve prediction IDs, ``agn_locus_refr_mrna_ids(locus)`` to retrieve reference IDs, or ``agn_locus_get_mrna_ids(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: GtUword agn_locus_mrna_num(AgnLocus *locus, AgnComparisonSource src)

  Get the number of mRNAs for the locus. Rather than calling this function directly, users are encouraged to use one of the following macros: ``agn_locus_num_pred_mrnas(locus)`` for the number of prediction mRNAs, ``agn_locus_num_refr_mrnas(locus)`` for the number of reference mRNAs, or ``agn_locus_num_mrnas(locus)`` if the source of annotation is undesignated or irrelevant.

.. c:function:: AgnLocus* agn_locus_new(GtStr *seqid)

  Class constructor.

.. c:function:: GtArray *agn_locus_pairs_to_report(AgnLocus *locus)

  Return the clique pairs to be reported for this locus.

.. c:function:: void agn_locus_png_track_selector(GtBlock *block, GtStr *track,void *data)

  Track selector function for generating PNG graphics of pairwise comparison loci. The track name to will be written to ``track``.

.. c:function:: void agn_locus_print_png(AgnLocus *locus, AgnLocusPngMetadata *metadata)

  Print a PNG graphic for this locus.

.. c:function:: void agn_locus_print_transcript_mapping(AgnLocus *locus, FILE *outstream)

  Print a mapping of the transcript(s) associated with this locus in a two-column tab-delimited format: ``transcriptId<tab>locusId``.

.. c:function:: void agn_locus_set_range(AgnLocus *locus, GtUword start, GtUword end)

  Set the start and end coordinates for this locus.

.. c:function:: double agn_locus_splice_complexity(AgnLocus *locus, AgnComparisonSource src)

  Calculate the splice complexity of this gene locus. Rather than calling this method directly, users are recommended to use one of the following macros: ``agn_locus_prep_splice_complexity(locus)`` to calculate the splice complexity of just the prediction transcripts, ``agn_locus_refr_splice_complexity(locus)`` to calculate the splice complexity of just the reference transcripts, and ``agn_locus_calc_splice_complexity(locus)`` to calculate the splice complexity taking into account all transcripts.

.. c:function:: bool agn_locus_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnLocusFilterStream
--------------------------

.. c:type:: AgnLocusFilterStream

  Implements the GenomeTools ``GtNodeStream`` interface. This is a node stream used to select loci based on user-specified criteria. See the `AgnLocusFilterStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocusFilterStream.h>`_.

.. c:function:: GtNodeStream* agn_locus_filter_stream_new(GtNodeStream *in_stream, GtArray *filters)

  Class constructor. The keys of the ``typestokeep`` hashmap should be the type(s) to be kept from the node stream. Any non-NULL value can be associated with those keys.

.. c:function:: bool agn_locus_filter_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnLocusMapVisitor
------------------------

.. c:type:: AgnLocusMapVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for printing out gene --> locus and mRNA --> locus relationships as part of a locus/iLocus processing stream. See the `AgnLocusMapVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocusMapVisitor.h>`_.

.. c:function:: GtNodeStream* agn_locus_map_stream_new(GtNodeStream *in, FILE *genefh, FILE *mrnafh)

  Constructor for a node stream based on this node visitor. See :c:func:`agn_locus_map_visitor_new` for a description of the function arguments.

.. c:function:: GtNodeVisitor *agn_locus_map_visitor_new(FILE *genefh, FILE *mrnafh)

  Constructor for the node visitor. Gene-to-locus relationships are printed to the ``genefh`` file handle, while mRNA-to-locus relationships are printed to the ``mrnafh`` file handle. Setting either file handle to NULL will disable printing the corresponding output.

Class AgnLocusStream
--------------------

.. c:type:: AgnLocusStream

  Implements the ``GtNodeStream`` interface. The only feature nodes delivered by this stream have type ``locus``, and the only direct children of these features are transcript features (of types mRNA, rRNA, or tRNA) present in the input stream. Any overlapping transcripts are children of the same locus feature. See the `AgnLocusStream class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocusStream.h>`_.

.. c:function:: GtNodeStream *agn_locus_stream_new(GtNodeStream *in_stream, GtLogger *logger)

  This constructor searches the complete feature graph of each feature node in the input stream for transcript features.

.. c:function:: GtNodeStream *agn_locus_stream_new_pairwise(GtNodeStream *refr_stream, GtNodeStream *pred_stream, GtLogger *logger)

  This constructor accepts two :c:type:`AgnTranscriptStream` objects as input. Locus features are created as per the class description, with additional data stored to track the source (reference vs prediction) of each transcript in each locus.

.. c:function:: void agn_locus_stream_set_idformat(AgnLocusStream *stream, const char *format)

  Loci created by this stream are assigned an ID with an arbitrary number. The default format is 'locus%lu' (that is, locus1, ;ocus2, etc). Use this function to override the default ID format.

.. c:function:: void agn_locus_stream_set_source(AgnLocusStream *stream, GtStr *source)

  Set the source value to be used for all iLoci created by this stream. Default value is 'AEGeAn::AgnLocusStream'.

.. c:function:: bool agn_locus_stream_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnMrnaRepVisitor
-----------------------

.. c:type:: AgnMrnaRepVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for filtering out all but the longest mRNA (as measured by CDS length) from alternatively spliced genes. See the `AgnMrnaRepVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnMrnaRepVisitor.h>`_.

.. c:function:: GtNodeStream* agn_mrna_rep_stream_new(GtNodeStream *in)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_mrna_rep_visitor_new()

  Constructor for the node visitor.

.. c:function:: bool agn_mrna_rep_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnNodeDeleteVisitor
--------------------------

.. c:type:: AgnNodeDeleteVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used to decrement the reference count to all feature nodes passing through the node stream. See the `AgnNodeDeleteVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnNodeDeleteVisitor.h>`_.

.. c:function:: GtNodeStream* agn_node_delete_stream_new(GtNodeStream *in)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_node_delete_visitor_new()

  Constructor for the node visitor.

Class AgnPseudogeneFixVisitor
-----------------------------

.. c:type:: AgnPseudogeneFixVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for correcting the ``type`` value for pseudogene features erroneously using the ``gene`` type instead of the more appropriate ``pseudogene`` type. See the `AgnPseudogeneFixVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnPseudogeneFixVisitor.h>`_.

.. c:function:: GtNodeStream* agn_pseudogene_fix_stream_new(GtNodeStream *in)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_pseudogene_fix_visitor_new()

  Constructor for the node visitor.

.. c:function:: bool agn_pseudogene_fix_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnRemoveChildrenVisitor
------------------------------

.. c:type:: AgnRemoveChildrenVisitor

  Implements the GenomeTools ``GtNodeVisitor`` interface. This is a node visitor used for correcting removing all children of each top-level feature. Psuedo-features are not modified. See the `AgnRemoveChildrenVisitor class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnRemoveChildrenVisitor.h>`_.

.. c:function:: GtNodeStream* agn_remove_children_stream_new(GtNodeStream *in)

  Constructor for a node stream based on this node visitor.

.. c:function:: GtNodeVisitor *agn_remove_children_visitor_new()

  Constructor for the node visitor.

.. c:function:: bool agn_remove_children_visitor_unit_test(AgnUnitTest *test)

  Run unit tests for this class. Returns true if all tests passed.

Class AgnTranscriptClique
-------------------------

.. c:type:: AgnTranscriptClique

  The purpose of the AgnTranscriptClique class is to store data pertaining to an individual maximal transcript clique. This clique may only contain a single transcript, or it may contain many. The only stipulation is that the transcripts do not overlap. Under the hood, each ``AgnTranscriptClique`` instance is a pseudo node (a GtFeatureNode object) with one or more transcript features as direct children. See the `AgnTranscriptClique class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTranscriptClique.h>`_.

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

.. c:function:: char *agn_transcript_clique_id(AgnTranscriptClique *clique)

  Retrieve the ID attribute of the transcript associated with this clique. User is responsible to free the string.

.. c:function:: GtArray *agn_transcript_clique_ids(AgnTranscriptClique *clique)

  Retrieve the ID attributes of all transcripts associated with this clique.

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

Functions for testing feature types. See the `AgnTypecheck module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTypecheck.h>`_.

.. c:function:: bool agn_typecheck_cds(GtFeatureNode *fn)

  Returns true if the given feature is a CDS; false otherwise.

.. c:function:: GtUword agn_typecheck_count(GtFeatureNode *fn, bool (*func)(GtFeatureNode *))

  Count the number of ``fn``'s children that have the given type.

.. c:function:: bool agn_typecheck_exon(GtFeatureNode *fn)

  Returns true if the given feature is an exon; false otherwise.

.. c:function:: bool agn_typecheck_gene(GtFeatureNode *fn)

  Returns true if the given feature is a gene; false otherwise.

.. c:function:: bool agn_typecheck_intron(GtFeatureNode *fn)

  Returns true if the given feature is an intron; false otherwise.

.. c:function:: bool agn_typecheck_mrna(GtFeatureNode *fn)

  Returns true if the given feature is an mRNA; false otherwise.

.. c:function:: bool agn_typecheck_pseudogene(GtFeatureNode *fn)

  Returns true if the given feature is declared as a pseudogene; false otherwise.

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

  Class used for unit testing of classes and modules. See the `AgnUnitTest class header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUnitTest.h>`_.

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

Collection of assorted functions that are otherwise unrelated. See the `AgnUtils module header <https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUtils.h>`_.

.. c:type:: AgnSequenceRegion

  This data structure combines sequence coordinates with a sequence ID to facilitate their usage together.



.. c:function:: GtArray* agn_array_copy(GtArray *source, size_t size)

  Similar to ``gt_array_copy``, except that array elements are treated as pointers and dereferenced before being added to the new array.

.. c:function:: double agn_calc_splice_complexity(GtArray *transcripts)

  Determine the splice complexity of the given set of transcripts.

.. c:function:: GtUword agn_feature_index_copy_regions(GtFeatureIndex *dest, GtFeatureIndex *src, bool use_orig, GtError *error)

  Copy the sequence regions from ``src`` to ``dest``. If ``use_orig`` is true, regions specified by input region nodes (such as those parsed from ``##sequence-region`` pragmas in GFF3) are used. Otherwise, regions inferred directly from the feature nodes are used.

.. c:function:: GtUword agn_feature_index_copy_regions_pairwise(GtFeatureIndex *dest, GtFeatureIndex *refrsrc, GtFeatureIndex *predsrc, bool use_orig, GtError *error)

  Copy the sequence regions from ``refrsrc`` and ``predsrc`` to ``dest``. If ``use_orig`` is true, regions specified by input region nodes (such as those parsed from ``##sequence-region`` pragmas in GFF3) are used. Otherwise, regions inferred directly from the feature nodes are used.

.. c:function:: void agn_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)

  Remove feature ``fn`` and all its subfeatures from ``root``. Analogous to ``gt_feature_node_remove_leaf`` with the difference that ``fn`` need not be a leaf feature.

.. c:function:: bool agn_feature_overlap_check(GtArray *feats)

  Returns true if any of the features in ``feats`` overlaps, false otherwise.

.. c:function:: GtUword agn_mrna_cds_length(GtFeatureNode *mrna)

  Determine the length of an mRNA's coding sequence.

.. c:function:: GtRange agn_multi_child_range(GtFeatureNode *top, GtFeatureNode *rep)

  If a top-level feature ``top`` contains a multifeature child (with multi representative ``rep``), use this function to get the complete range of the multifeature.

.. c:function:: int agn_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)

  Compare function for data type ``GtGenomeNode ``, needed for sorting ``GtGenomeNode `` stored in ``GtArray`` objects.

.. c:function:: int agn_sprintf_comma(GtUword n, char *buffer)

  Format the given non-negative number with commas as the thousands separator. The resulting string will be written to ``buffer``.

.. c:function:: int agn_string_compare(const void *p1, const void *p2)

  Dereference the given pointers and compare the resulting strings (a la ``strcmp``).

.. c:function:: GtStrArray* agn_str_array_union(GtStrArray *a1, GtStrArray *a2)

  Find the strings that are present in either (or both) of the string arrays.

