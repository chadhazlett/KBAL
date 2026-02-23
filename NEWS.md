# kbal Package News

## Version 0.1.1 (2024-09-10)
- Initial submission to CRAN.
- Updated dependencies to require R version 3.5.0 or higher.
- Fixed various documentation mismatches and notes.

## Version 0.1.2 (2025-03-12)
- Attempted fix for MKL-based SVD failures on CRAN checks.

## Version 0.1.3 (2025-07-05)
- Minor fix: fullSVD <- TRUE is set when falling back to full SVD after a truncated SVD failure.

## Version 0.1.4 (2026-02-23)
- Ensured that `base.weights` are consistently propagated through weight construction.
- Added `base.weights` support to `getw()`, `getdist()`, and `kbal()` to allow baseline sampling weights in balancing and discrepancy calculations.
