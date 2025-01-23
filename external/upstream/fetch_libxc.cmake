add_library(Libxc::libxc STATIC IMPORTED)

find_package(Libxc CONFIG QUIET
NO_CMAKE_PATH
NO_CMAKE_PACKAGE_REGISTRY
NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
)
if(TARGET Libxc::libxc)
  message(STATUS "Found LibXC: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable LibXC could not be located.")
endif()
