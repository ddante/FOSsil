cmake_minimum_required (VERSION 2.6)
project (FOS Fortran)

message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")

SET(EXECUTABLE_OUTPUT_PATH ../build)

SET(ELEMENT_FILES
  Src/Elements/Element_Class.f90
  Src/Elements/Face_Class.f90
  Src/Elements/Segment_Class.f90
  Src/Elements/Triangle_Class.f90
  Src/Elements/Quadrangle_Class.f90
)

SET(ALGEBRA_FILES
  Src/Algebra/Lin_Algebra.f90
)

SET(QUADRATURE_FILES
  Src/Quadrature/Quadrature_rules.f90
)


SET(SCHEME_FILES
  Src/Schemes/Num_Scheme.f90
  Src/Schemes/LDA_method.f90
  Src/Schemes/LW_method.f90
  Src/Schemes/LLxFS_method.f90
)

SET(MAIN_FILES 
  Src/Init_problem.f90
  Src/Models.f90
  Src/FOS_system.f90
  Src/Geometry.f90
  Src/Time_integration.f90
  Src/Space_integration.f90
  Src/Main.f90
  Src/Post_pro.f90
)

SET(SOURCES
  ${ELEMENT_FILES}
  ${ALGEBRA_FILES}
  ${QUADRATURE_FILES}
  ${SCHEME_FILES}
  ${MAIN_FILES}
)

ADD_EXECUTABLE(solver.exe ${SOURCES})

SET_TARGET_PROPERTIES(solver.exe  PROPERTIES COMPILE_FLAGS "-r8 -O5 -warn all")
#SET_TARGET_PROPERTIES(solver.exe  PROPERTIES COMPILE_FLAGS "-C -check uninit -ftrapuv -fpe0 -check noarg_temp_created -warn all")
