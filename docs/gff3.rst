Notes on GFF3
=============

AEGeAn's preferred exchange format
----------------------------------

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

*Ad hoc* conventions
--------------------

While the GFF3 specification requires the use of terms in the Sequence Ontology
for describing genomic features, it does not dictate the level of detail to
which genomic features must be described for a particular context. In many
cases, one could use different sets of terms to describe the same genomic
features. For example, features of a protein-coding gene include exons, introns,
a coding sequence, untranslated regions, a start codon, and a stop codon.
However, it is rare for all of these features to be declared explicitly in a
GFF3 file. Intron features are rarely included since they can easily be inferred
from exon structure. The coding sequence is sometimes encoded using a CDS
feature (defined with multiple entries if it spans multiple exons), but other
times is encoded using start and stop codons. Some GFF3 files list each mRNA as
a subfeature of an explicitly declared gene feature, while other files do not
explicitly declare gene features. The AEGeAn tools and API are designed to be
flexible and should work correctly so long as the features of interest are
described in sufficient detail.

Common pitfalls
---------------

* Some programs produce GFF3 that does not explicitly declare gene features for
  every transcript. The AEGeAn Toolkit includes code that can correct this issue
  by creating a gene parent with the same genomic location and attaching it as a
  parent to the transcript. However, this becomes problematic if, for example,
  the annotation includes alternatively spliced mRNAs without explicitly
  declared gene parents. Instead of inferring a single gene parent for all the
  alternative isoforms, AEGeAn will create a distinct gene feature for each
  individual isoform, leading to errors in subsequent analysis and reporting.
  If possible, users are much safer using GFF3 with gene features explicitly
  declared.
* Many whole-genome annotations, including those produced by NCBI's GNOMON
  annotation pipeline, include an annotation of organellar genomes such as
  those in mitochondria or chloroplasts. Users are responsible to decide which
  annotations to include in a particular analysis and should pre-filter their
  data accordingly. It is also important to note that genes from organellar
  genomes are sometimes encoded differently than genes from the nuclear genomes
  (this is true for GFF3 files produced by NCBI) and may need additional
  pre-processing before they can be analyzed with AEGeAn.
* The Sequence Ontology includes the terms `pseudogene` (to describe a
  pseudogene's location), `pseudogenic_exon` (to describe the pseudogene's
  structure), and `pseudogenic_transcript` in cases where the pseudogene is
  transcribed. However, some GFF3 files do not correctly use these terms (I'm
  looking at you NCBI) and use terms like `gene` and `exon` instead, leading to
  confusion and issues with downstream analysis. Fortunately, feature attributes
  can often be used to correct this issue. For example, pseudogenes in GFF3
  files produced by the `annotwriter` program from `NCBI's C++ Toolkit`_ label
  pseudogenes as genes, but include a `pseudo=true` attribute in the 9th column
  of the record. AEGeAn can correct this issue, but sometimes the user must
  request this behavior.

.. _`NCBI's C++ Toolkit`: http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC

References
----------

**Gremme G, Steinbiss S, Kurtz S** (2013) GenomeTools: A Comprehensive Software
Library for Efficient Processing of Structured Genome Annotations. *IEEE/ACM
Transactions on Computational Biology and Bioinformatics*, **10**:3, 645-656,
`doi:10.1109/TCBB.2013.68 <http://dx.doi.org/10.1109/TCBB.2013.68>`_.