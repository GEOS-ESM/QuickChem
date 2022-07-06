esma_set_this (OVERRIDE QC_Shared)

set (srcs
  QuickChem_Generic.F90
  xgb_fortran_api.F90
  )

find_package(xgboost 1.6.0 EXACT REQUIRED)

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL xgboost::xgboost)
