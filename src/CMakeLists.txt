# <<< Build MRChem library >>>
add_library(mrchem
  driver.cpp
  parallel.cpp
  mrenv.cpp
  )

target_include_directories(mrchem
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )

list(APPEND _public_headers
  driver.h
  mrchem.h
  mrenv.h
  parallel.h
  ${PROJECT_BINARY_DIR}/config.h
  )

add_subdirectory(chemistry)
add_subdirectory(analyticfunctions)
add_subdirectory(qmfunctions)
add_subdirectory(qmoperators)
add_subdirectory(mrdft)
add_subdirectory(scf_solver)
add_subdirectory(utils)

target_link_libraries(mrchem
  PUBLIC
    Eigen3::Eigen
    MRCPP::mrcpp
    XCFun::xcfun
    nlohmann_json::nlohmann_json
  )

set_target_properties(mrchem
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
    MACOSX_RPATH ON
    OUTPUT_NAME "mrchem"
    PUBLIC_HEADER "${_public_headers}"
  )

install(
  TARGETS
    mrchem
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  PUBLIC_HEADER
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MRChem
    COMPONENT dev
  )

# <<< Build MRChem executables >>>
add_executable(mrchem.x mrchem.cpp)

target_sources(mrchem.x
  PUBLIC
    $<BUILD_INTERFACE:$<JOIN:${CMAKE_CURRENT_LIST_DIR}/,${_public_headers}>>
  )

set_target_properties(mrchem.x
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${MRChem_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
  )

target_link_libraries(mrchem.x
  PRIVATE
    mrchem
    MRCPP::mrcpp
    XCFun::xcfun
    nlohmann_json::nlohmann_json
  )

# RPATH fixing
file(RELATIVE_PATH _rel ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR} ${CMAKE_INSTALL_PREFIX})
if(APPLE)
  set(_rpath "@loader_path/${_rel}")
else()
  set(_rpath "\$ORIGIN/${_rel}")
endif()
file(TO_NATIVE_PATH "${_rpath}/${CMAKE_INSTALL_LIBDIR}" MRChem_RPATH)

set_target_properties(mrchem.x
  PROPERTIES
    MACOSX_RPATH ON
    SKIP_BUILD_RPATH OFF
    BUILD_WITH_INSTALL_RPATH OFF
    INSTALL_RPATH "${MRChem_RPATH}"
    INSTALL_RPATH_USE_LINK_PATH ON
  )

install(
  TARGETS
    mrchem.x
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  )

add_subdirectory(initial_guess)
add_subdirectory(input)
