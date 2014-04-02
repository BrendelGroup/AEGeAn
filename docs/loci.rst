Loci as a Coordinate System
===========================

The concept of an *interval locus (iLocus)* was formulated as an organizational principle intended to facilitate reproducibility during the early stages of genome assembly and annotation.
Genome sequencing is now more affordable than ever, yet "finished" genome assemblies and annotations still require immense time and effort, and have become the exception more than the rule.
Embracing the new reality in which most genome projects will never progress beyond the preliminary draft stage, iLoci provide a robust coordinate system for working with rapidly changing genome annotation data.

.. image:: bogus.png
   :width: 500
   :alt: Awesome graphic here

Operational definition
----------------------

iLoci are intended to represent distinct, independent regions of a genome the encode one or more genes, as well as gene-less intergenic regions.
Given one or more pre-computed sets of gene annotations and parameter δ, iLoci are computed as follows.

* Create a bin for each gene in the input.
* If any of the genes in a bin overlap with any of the genes in another bin, merge those bins. Repeat until all possible merges are done.
* Sort the bins according to genomic position, and then determine the distance between each pair of adjacent bins.
    * If the distance is greater than 3δ, extend the bins toward each other by δ nucleotides.
    * If the distance is less than 3δ but greater than δ, extend the bins toward each other until they meet.
    * If the distance is less than δ, extend the bins toward each other as much as possible without overlapping with each other's genes (the extensions will overlap).
* Each bin now corresponds to a single iLocus. If there is empty space between any pair of adjacent iLoci, this also corresponds to an iLocus.

Note that while we define iLoci in terms of genes, the definition extends easily to any annotated genomic feature, and AEGeAn supports calculating iLoci for arbitrary feature types.
It's also important to keep in mind that iLoci are only as reliable as the annotations upon which they are based, and are subject to any technical artifacts present in those annotations.

Origins
-------

iLoci are based on the concept of a *gene locus*, which was originally formulated during development of the ParsEval program.
Given two sources of annotation for the same genomic sequence(s), the motivation behind the gene locus was to identify distinct protein-coding regions of the genome in which the two annotation sets could be analyzed and compared independently of all other regions.
ParsEval defines a gene locus as the smallest region containing all gene annotations that overlap with any other gene annotations in that region.
This definition ensured that no distinct gene loci overlap, that any gene annotation belongs only to a single gene locus, and that flanking intergenic regions can be ignored when computing statistical measures of similarity between the two sets of annotation at a particular gene locus.

The concept of an *iLocus* is simply an extension of this definition, with two main differences: the ends of each locus are extended using the δ parameter as described above, and intergenic regions are also considered iLoci.
