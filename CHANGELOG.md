# Change Log
All notable changes to this project will be documented in this CHANGELOG.md file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased Changes
- None.

## [0.9.0] - 2024-10-21
### Added
- Initial package release.

## [0.9.1] - 2024-10-30
### Added
- Better error message for `read_ensembl_gtf()` function that now points user to vignette for help when trying to read unsuported file types.
- Added `kaleido` as a dependency to allow users to images as multiple formats without having to install new packages.
- Improved code documentation, corrected spelling errors on vignettes, added two new vignettes.

## [0.9.2] - 2024-11-04
### Added
- Increased efficiency of to_intron().


## [0.9.3] - 2024-11-07
### Added
- Nothing... Just needed a new release to sync with Zotero for administrative reasons

## [0.9.4] - 2024-11-12
### Added
- Changed intron coordinates to be consistent with GTF format (1-index based, inclusive on both sides)
- Changed make traces function to make sure introns and exon/CDS "touch" in the plot even though their GTF coordinates do not touch.

## [0.9.5] - 2024-11-19
### Removed
- Kaleido dependency. Kaleido now requires chrome/chromium to be installed which could complicate the installation.


## [1.0.1] - 2025-02-10
### Added
- Changed make_traces() function to avoid transparent dots on the legend.
- Improved documentation for read_ensembl_gtf() function.
- Improved documentation for read_expression_matrix() function.
- Improved formatting on the vignettes.
- Improved README.md file.


## [1.0.2] - 2025-02-27
### Added
- Increase infer_schema_length when reading counts matrix to avoid type errors.



## [1.1.0] - 2025-03-12
### Added
- to_intron now automatically calculates exon number if it is not already present.
- shorten_gaps now automatically calculates exon_number if it is not already present.
- Adressed issue reported on https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/issues/6#issuecomment-2719061461, improving generalizability of the shorten gaps function
- Improved documentation for read_ensembl_gtf()

## [1.2.7] - 2025-03-27
### Added
- Correction to make the legend for expression match the order of the expression boxes.


## [1.2.8] - 2025-04-10
### Added
- Correction to the annotation legend to make it match the annotation order.


## [1.2.9] - 2025-04-23
### Added
- Removed debugging print statements

## [1.2.10] - 2025-04-24
### Added
- Further removed debugging print statements

## [1.2.11] - 2025-09-02
### Added
- Updated dependencies


## [1.2.12] - 2025-09-02
### Added
- Fixed tests