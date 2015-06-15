# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [0.13.1] - 2015-06-15

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
