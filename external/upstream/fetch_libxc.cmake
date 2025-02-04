cmake_policy(SET CMP0144 NEW)

add_library(Libxc::xc STATIC IMPORTED)
set_target_properties(Libxc::xc PROPERTIES
    IMPORTED_LOCATION "/home/ylva/work/libxc-mrchem/libxc/libxc_install/lib/libxc.a"
)


# set(export LIBXC_DIR=/home/ylva/work/libxc-mrchem/libxc/libxc_install)

find_package(Libxc REQUIRED CONFIG QUIET
NO_CMAKE_PATH
NO_CMAKE_PACKAGE_REGISTRY
NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
)
if(TARGET Libxc::xc)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Found LibXC: ${_loc} (found version ${Libxc_VERSION})")
else()
message(STATUS "Suitable LibXC could not be located. Fetching and building!")
include(FetchContent)
FetchContent_Declare(libxc_sources
  QUIET
  URL
    https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.gz
  )

FetchContent_GetProperties(libxc_sources)

set(CMAKE_BUILD_TYPE Release)
set(ENABLE_TESTALL FALSE CACHE BOOL "")
set(LIBXC_MAX_ORDER 3)  # TODO Maybe as a user-facing option?
set(LIBXC_PYTHON_INTERFACE FALSE CACHE BOOL "")

  if(NOT libxc_sources_POPULATED)
    FetchContent_Populate(libxc_sources)

    add_subdirectory(
      ${libxc_sources_SOURCE_DIR}
      ${libxc_sources_BINARY_DIR}
      )
  endif()
endif()
