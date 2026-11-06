# QuickChem

Quick Chemistry gridded component: a composite holding fast-chemistry species; today one child, hydroxyl radical (OH) predicted by an XGBoost model.

**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `QuickChem` @ `6547c1b` (2026-11-06).

**Build status on this line:** not built, not ported. On `release/MAPL-v3`, `@GEOSchem_GridComp/CMakeLists.txt` calls `esma_add_subdirectories(GOCART)` only, with the comment "Restore subdirs as components are ported to MAPL3."; `QuickChem` is not named in that file at all, not even in the commented-out list of subdirectories (`Shared GOCART TR GMI StratChem MAM MATRIX CARMA GAAS ACHEM`), and `GEOS_ChemGridComp.F90` has no reference to it. The code in this repository is the MAPL2 code (it differs from `main` only in `OH_GridComp/OH_instance_OH.rc`, the XGBoost model file path), kept for a later port. The Chemistry top-level library (`GEOS_ChemGridComp.F90`) is present but not compiled on this line: the `esma_add_library` block that builds it is commented out in `../CMakeLists.txt`. Statements about it describe the file as it is. This repository's own `CMakeLists.txt` lists the dependencies `MAPL`, `QC_Shared` (its own `Shared` subdirectory) and `esmf`; it lists neither `Chem_Base` nor `Chem_Shared`.

## Where it sits
GEOSgcm fixture (`components.yaml` @ `f50622d`, key `QuickChem`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@QuickChem`. Parent: Chemistry (`GEOSchem_GridComp`). **Not wired in**: the parent has no `MAPL_AddChild` for it and does not build it. Other fixtures may place it elsewhere.

## Documents
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): placement, file map, entry points, state, configuration, dependencies, traps, observations. Links up to the Chemistry architecture (`../docs/ARCHITECTURE.md` in the fixture).
- [docs/decisions.md](docs/decisions.md): architectural decisions.

## Layout
- `QuickChem_GridCompMod.F90`: composite (SetServices, Initialize, Run1, Run2), instance lists from `QuickChem_GridComp.rc`.
- `OH_GridComp/`: OH child (`OH_GridCompMod.F90`, `OH_StateSpecs.rc`, `OH_instance_OH.rc`, `OH_GridComp_ExtData.rc`).
- `QC_Environment/`: base type with `load_from_config`.
- `Shared/`: `QuickChem_Generic.F90`, `xgb_fortran_api.F90` (XGBoost 1.6.0 exactly).

## Rules for agents
- Note bugs and oddities; do not fix them unless asked.
- Keep the documents current when you change code; every factual claim carries a `(file, routine)` citation or an *Inferred:* / *Unverified:* label.
- Documentation changes only touch `AGENTS.md` and `docs/`; `CHANGELOG.md` is enforced by CI for code changes.
- Instance names containing `data` select the data-driven branch.
- `NOTES.wiki` is informal (2023) and unverified.

## Build and test
- Not built on this line (see build status above). Meant to be built as part of GEOSgcm by a parent that includes it (none does); needs XGBoost 1.6.0 and MAPL (`Shared/CMakeLists.txt`). No tests in this repository.
