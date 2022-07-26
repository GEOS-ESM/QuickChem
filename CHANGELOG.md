# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2022-07-26

### Changed
- Optimized OH by reducing the number of calls to XGBoost

### Added

- Initial release of QuickChem repository
- Added changelog enforcer and yaml validator GitHub actions
- Added routine IS_QC_INSTANCE_RUNNING, to be called by the parent GridComp
- Added OH flag compute_once_per_day, and support for RUNALARM
