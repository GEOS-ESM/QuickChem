# QuickChem decisions

**Audience:** coding agents and developers working on QuickChem.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `QuickChem` @ `6547c1b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).

### Instance-list composite copied from GOCART2G
- Pattern: the composite reads `ACTIVE_INSTANCES_<species>` / `PASSIVE_INSTANCES_<species>` from `QuickChem_GridComp.rc`, adds one child per name with `MAPL_AddChild`, active first (`QuickChem_GridCompMod.F90`, `getInstances_`, `createInstances_`). [Verified]
- Likely reason: reuse of the GOCART2G structure (revision history: "adapted from GOCART2G"); the ordering comment mentions "QC_AERO states that are passed to radiation" (`QuickChem_GridCompMod.F90`, `createInstances_`), and `QuickChem_GridComp.rc` says "Actually QuickChem does not provide aerosols!" (`AERO_PROVIDER` remark in its header). [Inferred]

### Data-or-computational by instance name
- Pattern: a child name containing `data` is a data-driven instance (`Shared/QuickChem_Generic.F90`, `determine_data_driven`); the composite skips Run2 for such children (`QuickChem_GridCompMod.F90`, Run2). [Verified]
- Likely reason: same convention as GOCART2G. [Inferred]

### OH predicted by XGBoost, once per day by default
- Pattern: OH comes from an XGBoost model (one file per month, `%m2` template) on a column of meteorology, composition and aerosol inputs; `compute_once_per_day: T` calls it only at 00:00:00; the previous result is kept in a private array (`OH_GridComp/OH_GridCompMod.F90`, Run1, `predict_OH_with_XGB`; `OH_GridComp/OH_instance_OH.rc`). [Verified]
- Likely reason: cost of the model; daily-mean inputs (`ONLINE_AVG24`) change slowly. [Inferred]

### Species work done in Run1, Run2 left as a stub
- Pattern: OH computes in phase 1; phase 2 only checks and resets the run alarm (`OH_GridCompMod.F90`, Run2). [Verified]
- Likely reason: placeholder structure from GOCART2G. [Inferred]
- Question for maintainers: is the reset in Run2 required by the alarm design (the Run1 comment says so)?

### Not wired into Chemistry
- Pattern: no `MAPL_AddChild` for QuickChem and no entry in the parent `CMakeLists.txt` (`../GEOS_ChemGridComp.F90`, `../CMakeLists.txt`); on `release/MAPL-v3` the parent builds only GOCART and does not name QuickChem even in its commented-out subdirectory list. [Verified]
- Question for maintainers: is QuickChem meant to become a Chemistry child, and when?
