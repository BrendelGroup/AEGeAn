# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]

### Changed
- Adopted strict adherance to the "Allman" C formatting style, as enforced by `clang-format`. (`clang-format -i -style=file src/pmrna.c`)
- Required C99 supped (only because C11 support isn't more ubiquitous).

## [0.15.1] - 2016-01-08

### Fixed
- Updated installation instructions.

## [0.15.0] - 2016-01-08

### Changed
- Now using a "label" (accession/Name/ID/position) in place of feature IDs in many places.
- iLocus children and grandchildren counts are now prefixed with `child_`, so as to prevent warning messages for features with an uppercase first letter
  (GFF3 attribute keys starting with an uppercase letter are disallowed except for those declared in the specification).
- Unit tests adjusted to handle GenomeTools 1.5.8's handling of the `gff-version` pragma (1 space instead of 3).

### Fixed
- Bug in `miloci.py` that prevented the last miLocus from printing.
- Better reporting of issues with gene/mRNA interval containment.
- Make build order deterministic (to facilitate reproducible builds on, i.e., Debian).
- Implemented Python 3 support for all ancillary Python scripts.
- A bug with the locus filtering code that distinguished between reference and prediction loci.
- Improved build and test environment for systems without Cairo installed.

## [0.14.1] - 2015-10-02

### Fixed
- Date in ChangeLog

## [0.14.0] - 2015-10-02

### Added
- New post-processing steps to iLocus parsing procedures (iLocus "refinement").
- An option for using coding sequence instead of UTRs as boundary test for gene overlap when computing iLoci.
- A variety of new iLocus GFF3 attributes related to iLocus accounting (such as for computing genome breakdown).
- A barrage of new functional tests (and some unit tests) to validate new functionality.
- Identification of intron genes into separate iLoci during iLocus parsing.

### Changed
- Replaced ParsEval's `--png` option with a `--nopng` option, made graphics the default for HTML output mode.
- iLocus delta extensions are now `\delta` bp long, even if that means extending into adjacent gene bodies.
- Protein-coding genes and non-coding genes are no longer placed in the same iLocus even if they overlap.
- Enabled use of feature accessions rather than IDs for some tools.
- Replaced the term *empty* with the more current term *iiLocus* in a variety of places in the code and documentation.
- Refactored iLocus functional test code.

### Fixed
- Bug that prevented AEGeAn from building correctly--thanks Sascha!
- Bug with option parsing and unknown options/flags (only appeared in Linux).

## [0.13.0] - 2015-06-15

### Added
- Ability to filter based on feature attribute key/value pairs.
- Ability to specify arbitrary parent types for mRNA rep stream, and corresponding `--locus` flag for the `pmrna` program.
- The `map` option for mRNA rep visitor and `pmrna`.
- Ability to access ParsEval locus reports by category (perfect match, CDS match, etc).
- Delta option for ParsEval.

### Fixed
- Simplified and corrected GAEVAL interval merging function.
- Fixed locus parsing issues resulting from incorrect sorting of input when pseudonodes involved.

## [0.12.2] - 2015-03-06

### Fixed
- Documentation issues (installation instructions and C API docs).

## [0.12.1] - 2015-03-06

### Fixed
- Change log formatting issues.
- Correct year in AgnVersion.h.

### Removed
- Faulty Coverity Scan configuration.

## [0.12.0] - 2015-03-05

### Fixed
- Fixed CDS phase correction on the reverse strand.
- Corrected non-deterministic selection of primary isoform when multiple isoforms have the same (maximum) length.
- In addition to correcting gene features, `AgnPseudogeneFixStream` also now corrects transcript and exon features as well.
- Corrected memory leak in pmrna.
- Implemented the advertised (but missing) `set_source` functionality in `AgnInferCDSVisitor` class.

### Changed
- Locations printed by Xtractore (of the form `seqid_start-end`) now include a strand character if strand is defined.
- Began using [this specification](http://keepachangelog.com/) for maintaining a ChangeLog.

### Added
- Python script for merging iLoci.
- Documentation for GFF3 expectations, both in general and tool-specific terms.
- Functional tests for primary mRNA selection and pseudogene correction.
- Node visitor for calculating coverage and GAEVAL integrity scores.
- New implementation of GAEVAL that calculates coverage and integrity scores for mRNA features.
- New branch `coverity`, which leverages Travis to run Coverity Scan static analysis.

### Removed
- Discarded vestigial `AgnNodeDeleteVisitor` class.
