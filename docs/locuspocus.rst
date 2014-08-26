LocusPocus
==========

Introduction
------------

**LocusPocus** is a program for computing interval loci (iLoci) from a provided
set gene annotations. iLoci correspond to a single gene, a set of overlapping
genes, or a space between genes. See :doc:`this page <loci>` for a description
of iLoci as an organizational principle for genomics.

Input
-----

Input for LocusPocus is one or more files in GFF3 format (LocusPocus can also
read from standard input if `-` is provided as the input filename). Aside from
compliance to `GFF3 syntax <http://sequenceontology.org/resources/gff3.html>`_,
the only requirement LocusPocus places on the input data is that gene features
must be explicitly declared. Note: LocusPocus output may be more descriptive
depending on whether additional features are described in the input, but the
iLoci themselves depend only on gene features.

Alternatively, users can override `gene` as the default feature of interest,
replace it with one or more feature types, and construct iLoci for these
features in the same way.

Output
------

LocusPocus computes the location of the iLoci from the given gene features and
reports the iLocus locations in GFF3 format. By default, only the locus feature
itself is reported, with attributes that indicate the number of genes in the
locus (as well as gene children such as mRNA, tRNA, etc, if they are provided
in the input). Invoking the `--verbose` option will additionally report all gene
features (and their subfeatures) associated with the locus.

Running LocusPocus
------------------

For a description of LocusPocus' command-line interface, run the following
command (after AEGeAn has been installed).

.. code-block:: bash

    locuspocus --help
