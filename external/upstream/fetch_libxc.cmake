cmake_policy(SET CMP0144 NEW)
find_package(Libxc CONFIG QUIET
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
  FetchContent_Declare(Libxc
    QUIET
    URL
      https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.gz
    )

  FetchContent_GetProperties(Libxc)

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_TESTALL FALSE CACHE BOOL "")
  set(LIBXC_MAX_ORDER 3)  # TODO Maybe as a user-facing option?
  set(LIBXC_PYTHON_INTERFACE FALSE CACHE BOOL "")

  if(NOT libxc_sources_POPULATED)
    FetchContent_Populate(Libxc)

    add_subdirectory(
      ${libxc_sources_SOURCE_DIR}
      ${libxc_sources_BINARY_DIR}
      )
  endif()
endif()
