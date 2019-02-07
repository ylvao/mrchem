if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
  add_executable(mrchem-pilot mrchem.cpp)
  target_include_directories(mrchem-pilot PRIVATE
      ${PROJECT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/src
      ${PROJECT_SOURCE_DIR}/src/chemistry
      ${PROJECT_SOURCE_DIR}/src/properties
      ${PROJECT_SOURCE_DIR}/src/analyticfunctions
      ${PROJECT_SOURCE_DIR}/src/qmfunctions
      ${PROJECT_SOURCE_DIR}/src/qmoperators
      ${PROJECT_SOURCE_DIR}/src/qmoperators/one_electron
      ${PROJECT_SOURCE_DIR}/src/qmoperators/two_electron
      ${PROJECT_SOURCE_DIR}/src/scf_solver
      )
  target_link_libraries(mrchem-pilot
      MRCPP::mrcpp
      XCFun::xcfun
      getkw::getkw
      )

  set(MRCHEM_DRIVER_DEBUG ON)
  set(MRCHEM_EXECUTABLE ${PROJECT_BINARY_DIR}/pilot/mrchem-pilot)

  configure_file(${PROJECT_SOURCE_DIR}/src/mrchem.in ${PROJECT_BINARY_DIR}/pilot/mrchem-pilot @ONLY)
  execute_process(COMMAND chmod 755 ${PROJECT_BINARY_DIR}/pilot/mrchem-pilot OUTPUT_QUIET)
endif()

if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.inp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
