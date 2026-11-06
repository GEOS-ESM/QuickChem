# QuickChem architecture

**Audience:** coding agents and developers working on QuickChem.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `QuickChem` @ e55c214 (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`, key `QuickChem`, local path `./src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSchem_GridComp/@QuickChem`, tag `v1.0.0`, develop branch `main`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths below are relative to this repository's root. Parent document: Chemistry architecture (`../../docs/ARCHITECTURE.md` in the fixture) (section 2.3 lists nested child repositories; section 3 describes the parent).

## 1. Placement

- Parent repository: `GEOSchem_GridComp` (the Chemistry component), at `../` from this repository. Its document records `@QuickChem` as "not wired in" (Chemistry architecture, section 2.3 (`../../docs/ARCHITECTURE.md` in the fixture)).
- Re-checked here: `grep -i quickchem` finds no match in `../GEOS_ChemGridComp.F90` (no `MAPL_AddChild` for it) and none in `../CMakeLists.txt`. So in the GEOSgcm fixture at `0f2e0fb`, QuickChem is checked out but neither compiled nor run by Chemistry. Nothing in this repository sets that up; the repository's own `CMakeLists.txt` is written to be included by a parent (`esma_set_this`, `esma_add_library`).
- *Inferred:* the module's header comment and `README.md` ("Quick Chemistry module") describe a component meant to sit among Chemistry's children; the exported field `OH` has a comment "Allow children of Chemistry to connect to these fields" (`QuickChem_GridCompMod.F90`, SetServices).
- No child repositories; no nested repositories.

## 2. Directory and file map

| Path | Content |
|---|---|
| `QuickChem_GridCompMod.F90` | module `QuickChem_GridCompMod`: the QuickChem composite gridded component; reads the instance list, adds the species children, runs them |
| `QuickChem_GridComp.rc` | instance lists: `ACTIVE_INSTANCES_OH`, `PASSIVE_INSTANCES_OH`, plus a placeholder species `XX` |
| `QuickChem_StateSpecs.rc` | automatic code generator (ACG) specifications for the composite itself: schema header only, no rows (empty IMPORT, EXPORT, INTERNAL tables) |
| `OH_GridComp/` | `OH_GridCompMod.F90` (module `OH_GridCompMod`, 1990 lines): the hydroxyl radical (OH) species child; `OH_StateSpecs.rc` (ACG specs); `OH_instance_OH.rc` (instance resources); `OH_GridComp_ExtData.rc` (ExtData rows, fragment; see 5) |
| `QC_Environment/` | `QC_EnvironmentMod.F90`: base type `QC_Environment` with `load_from_config` (reads `nbins:`); a `NOTE` file holds one question |
| `Shared/` | library `QC_Shared`: `QuickChem_Generic.F90` (`determine_data_driven`) and `xgb_fortran_api.F90` (Fortran interfaces to the XGBoost C library) |
| `NOTES.wiki` | informal setup notes dated March 2023 (see 9) |
| `CMakeLists.txt` (top level, `OH_GridComp/`, `QC_Environment/`, `Shared/`) | build; see 6 |
| `README.md`, `CHANGELOG.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `LICENSE`, `.github/` | standard repository files and continuous-integration workflows (`changelog-enforcer`, `enforce-labels`, `validate_yaml_files`) |

## 3. Purpose

QuickChem is a composite gridded component that holds "quick" chemistry species as children (`QuickChem_GridCompMod.F90`, module header: "includes OH, etc"). Only one species exists: OH (hydroxyl radical), either predicted online by a trained XGBoost model or read from a precomputed dataset (`OH_GridComp/OH_GridCompMod.F90`, `Run`, `Run1`, `Run_data`). It was adapted from GOCART2G (`QuickChem_GridCompMod.F90`, revision history "adapted from GOCART2G"); many GOCART2G remnants remain as comments.

## 4. Place in the hierarchy

- Parent: none today; Chemistry does not add it (section 1). `IS_QC_INSTANCE_RUNNING` (public) lets a parent ask whether a named instance is active, "to help parent GC determine whether to AddConnectivity" (`QuickChem_GridCompMod.F90`, `IS_QC_INSTANCE_RUNNING`). Nothing in the fixture calls it (grep in `@GEOSchem_GridComp`, excluding other nested repositories, finds no use).
- Children: one per name in `ACTIVE_INSTANCES_OH` then `PASSIVE_INSTANCES_OH`, each added with `MAPL_AddChild` and `OH_GridCompMod`'s `SetServices` (`QuickChem_GridCompMod.F90`, `getInstances_`, `createInstances_`, `addChildren__`). Default: one active instance `OH` (`QuickChem_GridComp.rc`). The `XX` lines in the rc file are read by nothing; the call for `XX` is commented out.
- Instance kind comes from the instance name: a name containing `data` is a data-driven instance, otherwise computational (`Shared/QuickChem_Generic.F90`, `determine_data_driven`), the same convention as GOCART2G (Chemistry architecture, section 3.1 (`../../docs/ARCHITECTURE.md` in the fixture)).
- The composite exports one field, `OH`, mapped to the first OH instance's export: `MAPL_AddExportSpec(GC, SHORT_NAME='OH', CHILD_ID=...)` when at least one instance exists (`QuickChem_GridCompMod.F90`, SetServices).

## 5. Entry points

Composite (`QuickChem_GridCompMod.F90`):
- SetServices: sets Initialize, Run phase 1 (`Run1`) and Run phase 2 (`Run2`) entry points; allocates a private state; reads `QuickChem_GridComp.rc`; adds children; includes `QuickChem_Export___.h`; adds the `OH` export; calls `MAPL_GenericSetServices`. No Finalize is set.
- Initialize: `MAPL_GenericInitialize`, fetches the private state and children; then loads `AGCM.rc` into a config object it does not use afterwards (comment: "should the rc file be changed?").
- `Run1`: runs every child, phase 1. `Run2`: runs every child whose name does not contain `data`, phase 2. Neither does work of its own (`Run2`'s comment records that GOCART2G's total diagnostics are not carried over).

OH child (`OH_GridComp/OH_GridCompMod.F90`):
- SetServices: loads `OH_instance_<name>.rc`, falling back to `OH_instance_OH.rc`; reads options (see 7); loads `GOCART2G_GridComp.rc` from the run directory to find the index of `wavelength_for_scacoef` in `wavelengths_for_profile_aop_in_nm:` (asserts if absent); sets Initialize, Run (phase 1) and, for computational instances only, Run2 (phase 2).
- Initialize: reads `RUN_DT:` and `QUICKCHEM_DT:` (default the model step), allocates `OH_ML`, calls `MAPL_GenericInitialize`.
- `Run` (phase 1): data instance calls `Run_data`; computational instance calls `Run1`.
- `Run1` (computational): returns at once unless the instance's run alarm is ringing. Otherwise assembles imports, calls `predict_OH_with_XGB` unless `compute_once_per_day` is true and the time of day is not 00:00:00, applies `OHscale`, and writes `OH` in the internal state: the XGBoost value where the layer pressure exceeds the tropopause pressure, `default_OH` elsewhere. Also fills the `DIAG_*` and `OH_boost` exports if requested.
- `Run2` (computational): checks the run alarm and switches it off (`ESMF_AlarmRingerOff`); no other work (GOCART2G processes are comments).
- `Run_data` (data instance): copies import `oh_OH` into the internal `OH`; asserts one bin.
- No Finalize (the entry point is commented out).

## 6. State

- Composite: ACG file `QuickChem_StateSpecs.rc` has no rows; the generated `QuickChem_Export___.h`, `QuickChem_DeclarePointer___.h`, `QuickChem_GetPointer___.h` are therefore empty of fields (*Inferred:* from the empty tables). The only export is `OH`, added by hand.
- OH child ACG file `OH_StateSpecs.rc`: internal `OH` (molec cm-3, xyz, restart optional, add2export true); imports `oh_*` (MERRA2-GMI "precomputed" fields: trace gases, `oh_ALBUV`, `oh_GMITO3`, `oh_OH`); exports `DIAG_*` diagnostics and `OH_boost` (mol/mol). Imports for the online sources (`T`, `Q`, `PLE`, `ZLE`, `TAUCLW`, `TAUCLI`, `TROPP`, aerosol scattering coefficients `*SCACOEF`, with `_avg24` daily-mean variants) are added by hand in SetServices by macros, depending on `OH_data_source` (below). Daily-mean imports use a 24-hour refresh and averaging interval and are not restarted (`MAPL_RestartSkip`); when `OH_data_source == ONLINE_AVG24`, Run1 gets `T_avg24` and sets `use_inst_values = .TRUE.` if `ptr3d(1,1,1) == 0.0` (`OH_GridCompMod.F90`, Run1). Its comment says that in the first 24 hours with no field in the import restart "the Coupler returns a field of all zero", so instantaneous imports are used, and that valid `T` is taken to mean all daily means are valid. The coupler behaviour is stated only in that comment; this repository does not show it (*Unverified:*).
- Internal `OH` has two declarations. Computational instances: generated from `OH_StateSpecs.rc` (`OH_Internal___.h`, included in SetServices only when not data-driven). Data instances: added by hand in SetServices (`MAPL_AddInternalSpec`, units `'???'`, `MAPL_RestartOptional`, `add2export=.true.`, `ungridded_dims=[nbins]`, i.e. a 4-D field), together with a `climoh001` import. The two declarations differ in units and in the ungridded dimension (*Inferred:* from the two code paths; `Run_data` copies a 3-D import into it).
- Key couplings (all through import specs; *Inferred:* the connections are made in the parent or in ExtData, which are not in this repository): `T`, `Q`, `PLE`, `ZLE`, `TAUCLW`, `TAUCLI`, `TROPP` from the model; `*SCACOEF` (4-D, wavelength dimension) from GOCART2G; `CH4`, `CO`, and the `oh_*` trace gases from ExtData or other children.
- `OH_GridComp_ExtData.rc` is an ExtData `PrimaryExports%%` fragment mapping `oh_*` names to MERRA2-GMI files under `/css/merra2gmi/pub/` plus an OMI reflectivity climatology; it also contains `OH_SRC` and `climxx001` rows copied from dust/aerosol (*Unverified:* use). It is installed under `etc/OH_GridComp`; this repository does not say who reads it.
- Restart: internal `OH` is optional-restart; the OH state kept between steps is the private array `OH_ML` (not saved in a restart) (`OH_GridCompMod.F90`, type `OH_GridComp`, Initialize).

## 7. Configuration that changes structure

- `QuickChem_GridComp.rc`: `ACTIVE_INSTANCES_OH`, `PASSIVE_INSTANCES_OH` decide how many OH children exist and their order; a name with `data` makes a data-driven instance.
- `OH_instance_OH.rc`: `nbins`; `XGBoostFile` (a grads template with `%m2`, one model file per month); `OH_data_source` (`PRECOMPUTED`, `ONLINE_INST`, `ONLINE_AVG24`; default file: `ONLINE_AVG24`) selects which imports exist and where inputs are read; `spinup_24hr_imports`; `compute_once_per_day`; `OHscale`; `wavelength_for_scacoef` (`OH_GridCompMod.F90`, SetServices).
- Run-time resources `QUICKCHEM_DT:` and `OH_DT`/`OH_REFERENCE_TIME` (the run alarm; comment block in Run1). *Unverified:* where `OH_DT` is read (by MAPL generic naming from the instance name).
- Compile time: XGBoost 1.6.0 exactly (`Shared/CMakeLists.txt`).

## 8. Dependencies

- MAPL, ESMF, NetCDF Fortran (`OH_GridComp/CMakeLists.txt`); XGBoost library `xgboost::xgboost`, version 1.6.0 exact, via `find_package` (`Shared/CMakeLists.txt`); the Fortran wrapper `Shared/xgb_fortran_api.F90` binds the C functions (`XGBoosterCreate`, `XGBoosterLoadModel`, `XGBoosterPredict`, `XGDMatrixCreateFromMat`, and others).
- Build structure: the top-level target is `QuickChem` with subcomponents `QC_Environment` and `OH_GridComp` and subdirectory `Shared` (library `QC_Shared`); `mapl_acg` runs on `QuickChem_StateSpecs.rc` and on `OH_StateSpecs.rc` (`CMakeLists.txt`, `OH_GridComp/CMakeLists.txt`). `OH_GridComp` installs every `*.rc` under it to `etc/OH_GridComp`; the top level installs `QuickChem_GridComp.rc` to `etc`.
- Run-time files: `QuickChem_GridComp.rc`, `OH_instance_<name>.rc`, `AGCM.rc`, `GOCART2G_GridComp.rc`, the XGBoost model files.
- Runtime cost: `predict_OH_with_XGB` builds a data matrix and loads the model each call (`xgb_fortran_api` calls; the free call is commented out, `OH_GridCompMod.F90`).

## 9. Known traps

- QuickChem is not wired into Chemistry (section 1); changes here have no effect on a GEOSgcm build until a parent adds it.
- OH requires `GOCART2G_GridComp.rc` in the run directory even when GOCART2G is not running (SetServices); the wavelength must be listed there.
- Data-driven instance: SetServices registers `oh_OH` and the other `oh_*` imports only inside the computational branch (the include `OH_Import___.h`); the data branch adds `climoh001` and `Run_data` reads `oh_OH`. See Observations.
- `QuickChem_GridComp.rc` header says "Default is the data component", but the default active instance is `OH` (computational).
- The instance name is a switch: any name containing the text `data` selects the data branch.

## 10. Decisions

See [decisions.md](decisions.md).

## 11. Observations (noted, not fixed) and open questions

- Data-driven path: `Run_data` reads the import `oh_OH`, whereas SetServices adds `climoh001` (not `oh_OH`) for data instances. *Inferred:* the `oh_*` import specs are added only for computational instances, so `MAPL_GetPointer(import, 'oh_OH')` may fail for a data instance. The comment "Is OH data driven? -- not implemented yet" is in SetServices. Question for maintainers: is the data-driven OH instance used?
- `OH_GridComp_ExtData.rc` carries `OH_SRC` (dust source) and `climxx001` rows that look inherited from GOCART2G; the OH child has no import with those names. Question: are they needed?
- QuickChem `Initialize` loads `AGCM.rc` into `CF` and does not use it again.
- `QuickChem_StateSpecs.rc` is empty yet `mapl_acg` runs on it; the generated headers are included.
- `Run2` of the composite and of OH do no chemistry; `OH` is computed in `Run1`, and `Run2` only silences the alarm. Question: is Run2 kept only as a placeholder?
- `XX` species placeholders remain in `QuickChem_GridComp.rc` and in the code as comments.
- `QC_Environment/NOTE` asks "Do we need a general routine like this to read the RC files?"
- Unit strings `???` in `OH_StateSpecs.rc` for most fields; the import `T`, `Q` rows in `OH_GridComp_ExtData.rc` list `Q` units as `Pa` (`oh_Q`).
- `NOTES.wiki` (Unverified, dated March 2023): describes building a feature branch of GEOSgcm and editing run files for an "ECCOH" experiment (adding CH4 to `Chem_Registry.rc`). Not checked against current code; not relied on above.
- `OH_REFERENCE_TIME` and `OH_DT` are described in comments only (`OH_GridCompMod.F90`, Run1).
