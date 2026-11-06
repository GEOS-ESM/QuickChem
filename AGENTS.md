# QuickChem

Quick Chemistry gridded component (MAPL2 line): a composite holding fast-chemistry species; today one child, hydroxyl radical (OH) predicted by an XGBoost model. Branch `main`, verified @ e55c214 (2026-11-06).

## Where it sits
GEOSgcm fixture (`components.yaml` @ `0f2e0fb`, key `QuickChem`): `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@QuickChem`. Parent: Chemistry (`GEOSchem_GridComp`). **Not wired in**: the parent has no `MAPL_AddChild` for it and does not build it. Other fixtures may place it elsewhere.

## Documents
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): placement, file map, entry points, state, configuration, dependencies, traps, observations. Links up to the [Chemistry architecture](../docs/ARCHITECTURE.md).
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
- Built only as part of GEOSgcm by a parent that includes it (it is not today); needs XGBoost 1.6.0 and MAPL (`Shared/CMakeLists.txt`). No tests in this repository.
- Use the NAG compiler first, then gfortran, via the `with-nag` and `with-gfortran` wrappers (workspace rules).
