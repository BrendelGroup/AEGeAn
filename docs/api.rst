AEGeAn C API
============

The most comprehensive documentation for AEGeAn's C API resides in the
project's header files. This page simply lists a brief description for each
class/module, as well as the methods associated with each. Please refer to the
header files (links provided) for more information on individual methods and
their arguments.

Class AgnCanonGeneStream
------------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCanonGeneStream.h.

A node stream that returns only canonical protein-coding genes that pass stringent validation.

Class AgnCliquePair
-------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnCliquePair.h.

The purpose of the AgnCliquePair class is to associate all of the data needed to compare transcripts from two alternative sources of annotation for the same sequence.

* ``AgnCliquePairClassification agn_clique_pair_classify(AgnCliquePair *pair)``

* ``void agn_clique_pair_comparative_analysis(AgnCliquePair *pair)``

* ``int agn_clique_pair_compare(void *p1, void *p2)``

* ``int agn_clique_pair_compare_direct(AgnCliquePair *p1, AgnCliquePair *p2)``

* ``int agn_clique_pair_compare_reverse(void *p1, void *p2)``

* ``void agn_clique_pair_delete(AgnCliquePair *pair)``

* ``double agn_clique_pair_get_edit_distance(AgnCliquePair *pair)``

* ``AgnTranscriptClique *agn_clique_pair_get_pred_clique(AgnCliquePair *pair)``

* ``const char *agn_clique_pair_get_pred_vector(AgnCliquePair *pair)``

* ``AgnTranscriptClique *agn_clique_pair_get_refr_clique(AgnCliquePair *pair)``

* ``const char *agn_clique_pair_get_refr_vector(AgnCliquePair *pair)``

* ``AgnComparison *agn_clique_pair_get_stats(AgnCliquePair *pair)``

* ``bool agn_clique_pair_has_utrs(AgnCliquePair *pair)``

* ``bool agn_clique_pair_is_simple(AgnCliquePair *pair)``

* ``unsigned long agn_clique_pair_length(AgnCliquePair *pair)``

* ``bool agn_clique_pair_needs_comparison(AgnCliquePair *pair)``

* ``AgnCliquePair* agn_clique_pair_new(const char *seqid, AgnTranscriptClique *refr_clique, AgnTranscriptClique *pred_clique, GtRange *locus_range)``

* ``void agn_clique_pair_record_characteristics(AgnCliquePair *pair, AgnCompResultDesc *desc)``

* ``bool agn_clique_pair_unit_test(AgnUnitTest *test)``

Module AgnComparEval
--------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnComparEval.h.

FIXME

* ``void agn_comp_evaluation_init(AgnCompEvaluation *data)``

* ``void agn_comp_result_summary_combine(AgnCompResultSummary *desc, AgnCompResultSummary *desc_to_add)``

* ``void agn_comp_result_summary_init(AgnCompResultSummary *desc)``

* ``void agn_comp_result_desc_combine(AgnCompResultDesc *desc, AgnCompResultDesc *desc_to_add)``

* ``void agn_comp_result_desc_init(AgnCompResultDesc *desc)``

* ``void agn_comp_summary_combine(AgnCompSummary *s1, AgnCompSummary *s2)``

* ``void agn_comp_summary_init(AgnCompSummary *summary)``

* ``void agn_comparison_combine(AgnComparison *c1, AgnComparison *c2)``

* ``void agn_comparison_init(AgnComparison *comparison)``

* ``void agn_compare_filters_init(AgnCompareFilters *filters)``

* ``void agn_compare_filters_parse(AgnCompareFilters *filters, FILE *instream, AgnLogger *logger)``

* ``void agn_comp_stats_binary_init(AgnCompStatsBinary *stats)``

* ``void agn_comp_stats_binary_resolve(AgnCompStatsBinary *stats)``

* ``void agn_comp_stats_scaled_init(AgnCompStatsScaled *stats)``

* ``void agn_comp_stats_scaled_resolve(AgnCompStatsScaled *stats)``

Class AgnFilterStream
---------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnFilterStream.h.

A node stream used to filter out unwanted nodes. The user can either specify the types of nodes to keep (all others will be deleted and not passed), or the types of nodes not to keep (these will be deleted, all others will be passed).

Class AgnGeneLocus
------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnGeneLocus.h.

The purpose of the AgnGeneLocus class is to store all of the data associated with a distinct locus, in many cases to facilitate the comparison of two sets of gene structure annotations for that locus.

* ``void agn_gene_locus_aggregate_results(AgnGeneLocus *locus, AgnCompEvaluation *eval)``

* ``AgnGeneLocus *agn_gene_locus_clone(AgnGeneLocus *locus)``

* ``int agn_gene_locus_array_compare(const void *p1, const void *p2)``

* ``unsigned long agn_gene_locus_cds_length(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``GtArray *agn_gene_locus_comparative_analysis(AgnGeneLocus *locus)``

* ``void agn_gene_locus_delete(AgnGeneLocus *locus)``

* ``unsigned long agn_gene_locus_exon_num(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``bool agn_gene_locus_filter(AgnGeneLocus *locus, AgnCompareFilters *filters)``

* ``GtArray *agn_gene_locus_genes(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``GtArray *agn_gene_locus_gene_ids(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``unsigned long agn_gene_locus_gene_num(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``unsigned long agn_gene_locus_get_end(AgnGeneLocus *locus)``

* ``unsigned long agn_gene_locus_get_length(AgnGeneLocus *locus)``

* ``const char* agn_gene_locus_get_seqid(AgnGeneLocus *locus)``

* ``unsigned long agn_gene_locus_get_start(AgnGeneLocus *locus)``

* ``GtArray *agn_gene_locus_get_unique_pred_cliques(AgnGeneLocus *locus)``

* ``GtArray *agn_gene_locus_get_unique_refr_cliques(AgnGeneLocus *locus)``

* ``AgnGeneLocus* agn_gene_locus_new(const char *seqid)``

* ``unsigned long agn_gene_locus_num_clique_pairs(AgnGeneLocus *locus)``

* ``void agn_gene_locus_png_track_selector(GtBlock *block, GtStr *track,void *data)``

* ``void agn_gene_locus_print_png(AgnGeneLocus *locus, AgnGeneLocusPngMetadata *metadata)``

* ``GtRange agn_gene_locus_range(AgnGeneLocus *locus)``

* ``void agn_gene_locus_set_range(AgnGeneLocus *locus, unsigned long start, unsigned long end)``

* ``double agn_gene_locus_splice_complexity(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``void agn_gene_locus_summary_init(AgnGeneLocusSummary *summary)``

* ``void agn_gene_locus_to_gff3(AgnGeneLocus *locus, FILE *outstream, const char *source)``

* ``GtArray *agn_gene_locus_transcripts(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``GtArray *agn_gene_locus_transcript_ids(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``unsigned long agn_gene_locus_transcript_num(AgnGeneLocus *locus, AgnComparisonSource src)``

* ``bool agn_gene_locus_unit_test(AgnUnitTest *test)``

Module AgnGtExtensions
----------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnGtExtensions.h.

A collection of extensions to core GenomeTools classes.

* ``void agn_gt_feature_index_to_gff3(GtFeatureIndex *index, FILE *outstream)``

* ``unsigned long agn_gt_feature_node_cds_length(GtFeatureNode *transcript)``

* ``GtArray *agn_gt_feature_node_children_of_type(GtFeatureNode *fn, bool (*typetestfunc)(GtFeatureNode *))``

* ``bool agn_gt_feature_node_is_cds_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_exon_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_gene_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_intron_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_mrna_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_start_codon_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_stop_codon_feature(GtFeatureNode *fn)``

* ``bool agn_gt_feature_node_is_utr_feature(GtFeatureNode *fn)``

* ``unsigned long agn_gt_feature_node_num_transcripts(GtFeatureNode *gene)``

* ``bool agn_gt_feature_node_overlap(GtFeatureNode *first, GtFeatureNode *second)``

* ``bool agn_gt_feature_node_range_contains(GtFeatureNode *n1, GtFeatureNode *n2)``

* ``void agn_gt_feature_node_remove_tree(GtFeatureNode *root, GtFeatureNode *fn)``

* ``void agn_gt_feature_node_resolve_pseudo_node(GtFeatureNode *root, GtArray *nodes)``

* ``void agn_gt_feature_node_set_source_recursive( GtFeatureNode *feature, GtStr *source )``

* ``void agn_gt_feature_node_to_gff3(GtFeatureNode *feature, FILE *outstream, bool printchildren, char *prefix, GtHashmap *filtered_types)``

* ``int agn_gt_genome_node_compare(const void *n1, const void *n2)``

* ``char agn_gt_phase_to_char(GtPhase phase)``

* ``char agn_gt_strand_to_char(GtStrand strand)``

* ``GtStrArray* agn_gt_str_array_intersection(GtStrArray *a1, GtStrArray *a2)``

* ``GtStrArray* agn_gt_str_array_union(GtStrArray *a1, GtStrArray *a2)``

Class AgnInferCDSVisitor
------------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferCDSVisitor.h.

FIXME

Class AgnInferExonsVisitor
--------------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnInferExonsVisitor.h.

FIXME

Class AgnLocusIndex
-------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLocusIndex.h.

FIXME

* ``void agn_locus_index_comparative_analysis(AgnLocusIndex *idx, const char *seqid, int numprocs, AgnLocusIndexVisitFunc preanalyfunc, AgnLocusIndexVisitFunc postanalyfunc, void *analyfuncdata, AgnLogger *logger)``

* ``void agn_locus_index_delete(AgnLocusIndex *idx)``

* ``void agn_locus_index_find(AgnLocusIndex *idx, const char *seqid, GtRange *range, GtArray *loci)``

* ``GtArray *agn_locus_index_get(AgnLocusIndex *idx, const char *seqid)``

* ``GtArray *agn_locus_index_interval_loci(AgnLocusIndex *idx, const char *seqid, unsigned long delta, bool skipterminal)``

* ``AgnLocusIndex *agn_locus_index_new(bool freeondelete)``

* ``unsigned long agn_locus_index_parse_pairwise_memory(AgnLocusIndex *idx, GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, int numprocs, AgnCompareFilters *filters, AgnLogger *logger)``

* ``unsigned long agn_locus_index_parse_pairwise_disk(AgnLocusIndex *idx, const char *refrfile, const char *predfile, int numprocs, AgnCompareFilters *filters, AgnLogger *logger)``

* ``unsigned long agn_locus_index_parse_memory(AgnLocusIndex *idx, GtFeatureIndex *features, int numprocs, AgnLogger *logger)``

* ``unsigned long agn_locus_index_parse_disk(AgnLocusIndex *idx, int numfiles, const char **filenames, int numprocs, AgnLogger *logger)``

* ``GtStrArray *agn_locus_index_seqids(AgnLocusIndex *idx)``

Class AgnLogger
---------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnLogger.h.

The AgnLogger class is desiged to store error, warning, and status messages.

* ``GtArray *agn_logger_get_error_messages(AgnLogger *logger)``

* ``GtArray *agn_logger_get_status_messages(AgnLogger *logger)``

* ``GtArray *agn_logger_get_warning_messages(AgnLogger *logger)``

* ``bool agn_logger_has_error(AgnLogger *logger)``

* ``bool agn_logger_has_status(AgnLogger *logger)``

* ``bool agn_logger_has_warning(AgnLogger *logger)``

* ``void agn_logger_log_error(AgnLogger *logger, const char *format, ...)``

* ``void agn_logger_log_status(AgnLogger *logger, const char *format, ...)``

* ``void agn_logger_log_warning(AgnLogger *logger, const char *format, ...)``

* ``AgnLogger *agn_logger_new()``

* ``bool agn_logger_print_all(AgnLogger *logger, FILE *outstream, const char *format, ...)``

* ``bool agn_logger_print_error(AgnLogger *logger, FILE *outstream, const char *format, ...)``

* ``bool agn_logger_print_status(AgnLogger *logger, FILE *outstream, const char *format, ...)``

* ``bool agn_logger_print_warning(AgnLogger *logger, FILE *outstream, const char *format, ...)``

* ``void agn_logger_unset(AgnLogger *logger)``

Class AgnTranscriptClique
-------------------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnTranscriptClique.h.

The purpose of the AgnTranscriptClique class is to store data pertaining to an individual maximal transcript clique. This clique may only contain a single transcript, or it may contain many. The only stipulation is that the transcripts do not overlap.

* ``void agn_transcript_clique_add(AgnTranscriptClique *clique, GtFeatureNode *transcript)``

* ``unsigned long agn_transcript_clique_cds_length(AgnTranscriptClique *clique)``

* ``AgnTranscriptClique* agn_transcript_clique_copy(AgnTranscriptClique *clique)``

* ``void agn_transcript_clique_delete(AgnTranscriptClique *clique)``

* ``bool agn_transcript_clique_has_id_in_hash(AgnTranscriptClique *clique, GtHashmap *map)``

* ``const char *agn_transcript_clique_id(AgnTranscriptClique *clique)``

* ``AgnTranscriptClique* agn_transcript_clique_new()``

* ``unsigned long agn_transcript_clique_num_exons(AgnTranscriptClique *clique)``

* ``unsigned long agn_transcript_clique_num_utrs(AgnTranscriptClique *clique)``

* ``void agn_transcript_clique_print_ids(AgnTranscriptClique *clique, FILE *outstream)``

* ``void agn_transcript_clique_put_ids_in_hash(AgnTranscriptClique *clique, GtHashmap *map)``

* ``unsigned long agn_transcript_clique_size(AgnTranscriptClique *clique)``

* ``GtArray* agn_transcript_clique_to_array(AgnTranscriptClique *clique)``

* ``void agn_transcript_clique_to_gff3(AgnTranscriptClique *clique, FILE *outstream, const char *prefix)``

* ``void agn_transcript_clique_traverse(AgnTranscriptClique *clique, AgnCliqueVisitFunc func, void *funcdata)``

* ``bool agn_transcript_clique_unit_test(AgnUnitTest *test)``

Class AgnUnitTest
-----------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUnitTest.h.

Class used for unit testing of classes and modules.

* ``void agn_unit_test_print(AgnUnitTest *test, FILE *outstream)``

* ``void agn_unit_test_result(AgnUnitTest *test, const char *label, bool success)``

* ``void agn_unit_test_run(AgnUnitTest *test)``

Module AgnUtils
---------------
See https://github.com/standage/AEGeAn/blob/master/inc/core/AgnUtils.h.

A collection of assorted utility functions for ParsEval.

* ``double agn_calc_splice_complexity(GtArray *transcripts)``

* ``GtFeatureNode *agn_eden()``

* ``GtArray* agn_enumerate_feature_cliques(GtArray *feature_set)``

* ``GtArray* agn_feature_neighbors(GtGenomeNode *feature, GtArray *feature_set)``

* ``FILE *agn_fopen(const char *filename, const char *mode, FILE *errstream)``

* ``GtFeatureIndex *agn_import_canonical(int numfiles, const char **filenames, AgnLogger *logger)``

* ``GtFeatureIndex *agn_import_simple(int numfiles, const char **filenames, char *type, AgnLogger *logger)``

* ``bool agn_infer_cds_range_from_exon_and_codons( GtRange *exon_range, GtRange *leftcodon_range, GtRange *rightcodon_range, GtRange *cds_range )``

* ``GtStrArray* agn_seq_intersection(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, AgnLogger *logger)``

* ``GtStrArray* agn_seq_union(GtFeatureIndex *refrfeats, GtFeatureIndex *predfeats, AgnLogger *logger)``

* ``int agn_sprintf_comma(unsigned long n, char *buffer)``

* ``int agn_string_compare(const void *p1, const void *p2)``

* ``GtRange agn_transcript_cds_range(GtFeatureNode *transcript)``

* ``void agn_transcript_structure_gbk(GtFeatureNode *transcript, FILE *outstream)``

