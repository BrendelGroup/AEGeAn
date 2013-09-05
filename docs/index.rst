.. AEGeAn documentation master file, created by
   sphinx-quickstart on Thu Sep  5 09:18:30 2013.

AEGeAn Toolkit: integrated genome analysis
==========================================

.. image:: aegean-logo-v1.0-128.png
   :align: right
   :alt: AEGeAn logo

:Author:  Daniel S. Standage
:Contact: daniel.standage@gmail.com
:License: ISC


**Contents:**

.. toctree::
   :maxdepth: 1
   
   install
   citing
   api
   license
   contrib


**Note**: if you're in a hurry, check out the :doc:`installation demo for the
impatient <install>`.

The AEGeAn Toolkit started as several distinct but related efforts to build
tools for managing and analyzing whole-genome gene structure annotations. AEGeAn
has brought these efforts together into a single library that includes
executable programs as well as several data structures and modules available via
a C API. The AEGeAn Toolkit leverages a variety of parsers, data structures, and
graphics capabilities available from the GenomeTools library
(http://genometools.org).

* **ParsEval** is a program for comparing distinct sets of gene structure
  annotations for the same sequence(s). This program calculates and reports a
  rich set of comparison statistics, both at the level of individual gene
  loci as well as at the level of entire sequences.

* **CanonGFF3** is a tool for preprocessing GFF3 data. It validates features
  related to canonical protein-coding genes, accepting data encoded in a wide
  variety of common conventions.

* **LocusPocus** is a program for computing gene loci from one or more gene
  prediction sets. In the ParsEval paper cited below, a 'gene locus' is
  defined as the smallest genomic region that contains all genes that overlap
  with any other genes in that region. This definition can be useful when
  comparing two sets of gene predictions.

* **VAnG** (Validation of Annotated Genomes) is a program for validating GFF3
  files against a schema.

* Additional tools are under development and will be released once they are a
  bit more stable.

If you have any questions regarding AEGeAn, feel free to contact the author by
email or, even better, open up a thread on `AEGeAn's issue tracker
<https://github.com/standage/AEGeAn/issues>`_ so that my response will be
visible to others who may have the same questions or issues in the future.
