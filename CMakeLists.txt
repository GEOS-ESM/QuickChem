###  based on @QuickChem/ESMF/QuickChem_GridComp/CMakeLists.txt :

esma_set_this ()

set (sources
  QuickChem_GridCompMod.F90
  )

set (resource_files
   QuickChem_GridComp.rc
   )

install( FILES ${resource_files}
   DESTINATION etc
   )

set (public_libs
  MAPL
  QC_Shared
  esmf
  )

esma_add_library (${this}
  SOURCES ${sources}
  SUBCOMPONENTS QC_Environment OH_GridComp
  SUBDIRS Shared
  PUBLIC_LIBS ${public_libs}
  )

mapl_acg (${this}   QuickChem_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)

