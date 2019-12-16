file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

configure_file (
  ${PROJECT_SOURCE_DIR}/config.h.in
  ${PROJECT_BINARY_DIR}/config.h
  @ONLY
  )

add_custom_command(
  OUTPUT
    ${PROJECT_BINARY_DIR}/version.h
  COMMAND
    ${CMAKE_COMMAND} -DINPUT_DIR=${PROJECT_SOURCE_DIR}
                     -DTARGET_DIR=${PROJECT_BINARY_DIR}
                     -DCMAKE_SYSTEM=${CMAKE_SYSTEM}
                     -DCMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}
                     -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
                     -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                     -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                     -DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                     -DVERSION_FILE=${PROJECT_SOURCE_DIR}/VERSION
                     -DMW_FILTER_SOURCE_DIR=${MW_FILTER_SOURCE_DIR}
                     -DMW_FILTER_INSTALL_DIR=${MW_FILTER_INSTALL_DIR}
                     -P ${CMAKE_CURRENT_LIST_DIR}/binary-info.cmake
  MAIN_DEPENDENCY
    ${PROJECT_SOURCE_DIR}/version.h.in
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_LIST_DIR}
  )

# rebuild version_info.h every time
add_custom_target(
  mrchem-info
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E touch_nocreate ${PROJECT_SOURCE_DIR}/version.h.in
  DEPENDS
    ${PROJECT_BINARY_DIR}/version.h
  )

# See here for the reason why: https://gitlab.kitware.com/cmake/cmake/issues/18399
set_source_files_properties(${PROJECT_BINARY_DIR}/version.h
  PROPERTIES
    GENERATED 1
  )

set(PYTHON_SITE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

# Fetch dependencies
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_xcfun.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_mrcpp.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_nlohmann_json.cmake)

add_subdirectory(src)
add_subdirectory(pilot)
