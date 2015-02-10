Notes on GFF3
=============

The AEGeAn Toolkit, and the GenomeTools library upon which AEGeAn heavily
relies, use the `GFF3 format`_ as the preferred encoding for genome annotations.
GFF3 utilizes a tab-delimited plain text format that bears resemblance, and
indeed owes its origins, to several related annotation encoding schemes. The
GFF3 format is unique, however, in at least two critical ways. First, it
leverages a controlled vocabulary, the `Sequence Ontology`_, for describing
genomic features and the relationships between them. Second, GFF3 enables
greater flexibility and granularity in describing genomic features compared to
alternative formats, such as the `Gene Transfer Format (GTF)`_ which focuses
exclusively on protein-coding genes, or the `Browser Extensible Data (BED)`_
format which focuses primarily on feature visualization.Both GTF and BED support
only a single level of feature decomposition, while GFF3 supports grouping
features and subfeatures to arbitrary levels of granularity. Annotations encoded
in GFF3 can be considered *annotation graphs*, directed acyclic graphs with
nodes representing genomic features and edges representing relationships between
the features (Gremme, 2013).

.. _`GFF3 format`: http://sequenceontology.org/resources/gff3.html
.. _`Sequence Ontology`: http://sequenceontology.org
.. _`Gene Transfer Format (GTF)`: http://mblab.wustl.edu/GTF22.html
.. _`Browser Extensible Data (BED)`: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

References
----------

**Gremme G, Steinbiss S, Kurtz S** (2013) GenomeTools: A Comprehensive Software
Library for Efficient Processing of Structured Genome Annotations. *IEEE/ACM
Transactions on Computational Biology and Bioinformatics*, **10**:3, 645-656,
`doi:10.1109/TCBB.2013.68 <http://dx.doi.org/10.1109/TCBB.2013.68>`_.