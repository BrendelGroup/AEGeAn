CanonGFF3
=========

Introduction
------------

**CanonGFF3** is a program for pre-processing GFF3 data encoding canonical
protein-coding genes. It will clean up a GFF3 file, removing all features not
directly related to protein-coding genes and inferring features that are not
explicitly declared, such as inrons and UTRs. Under the hood, CanonGFF3
essentially applies the same procedure used by ParsEval when it inspects its
GFF3 input.

Input
-----

Input for CanonGFF3 is one or more files in GFF3 format (CanonGFF3 can also read
from standard input if `-` is provided as the input filename). Aside from
compliance to `GFF3 syntax <http://sequenceontology.org/resources/gff3.html>`_,
CanonGFF3 requires only that protein coding genes be described in enough detail
that the entire gene structure can be interpreted. For example, one common
convention is to use exon and CDS features to describe the structure. No intron,
UTR, or start/stop features are *explicitly* provided, but these can be inferred
from the other features. An alternative convention is to only declare exon and
start & stop codon features, which requires the introns, UTRs, and CDS to be
inferred. CanonGFF3 is pretty flexible in its handling of these various
conventions, assuming the gene structure is described in sufficient detail.

Output
------

CanonGFF3 output is a GFF3 file containing protein-coding genes from the
provided input file(s). In most cases the output will be more verbose than the
input, containing features that have been inferred from the features provided
explicitly in the input.

Running CanonGFF3
-----------------

For a description of CanonGFF3's command-line interface, run the following
command (after AEGeAn has been installed).

.. code-block:: bash

    canon-gff3 --help
