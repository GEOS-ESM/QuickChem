esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES QC_Environment MAPL QC_Shared esmf NetCDF::NetCDF_Fortran)

mapl_acg (${this}   OH_StateSpecs.rc 
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS 
          GET_POINTERS DECLARE_POINTERS)

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
