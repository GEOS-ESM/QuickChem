# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added a scaling factor to scale OH by a constant
- Added exports for OH

### Changed

- Moved to GitHub Actions for label enforcement
- Changed OH Boost file to be climatology rather than single month

### Fixed

### Removed

## [1.0.0] - 2022-07-26

### Changed

- Optimized OH by reducing the number of calls to XGBoost

### Added

- Initial release of QuickChem repository
- Added changelog enforcer and yaml validator GitHub actions
- Added routine IS_QC_INSTANCE_RUNNING, to be called by the parent GridComp
- Added OH flag compute_once_per_day, and support for RUNALARM
