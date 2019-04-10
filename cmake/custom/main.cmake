configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

set(PYTHON_SITE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

# Fetch dependencies
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_xcfun.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_eigen3.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_mrcpp.cmake)
include(${PROJECT_SOURCE_DIR}/external/upstream/fetch_nlohmann_json.cmake)

add_subdirectory(src)
add_subdirectory(pilot)
