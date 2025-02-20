#set(Libxc_pinned "7.0.0")
cmake_policy(SET CMP0144 NEW)
find_package(Libxc CONFIG QUIET COMPONENTS C
NO_CMAKE_PATH
NO_CMAKE_PACKAGE_REGISTRY
NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
)
if(TARGET Libxc::xc)
  get_property(_loc TARGET Libxc::xc PROPERTY LOCATION)
  message(STATUS "Found LibXC: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable Libxc could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(Libxc
    QUIET
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    URL
      # https://gitlab.com/libxc/libxc/-/archive/${Libxc_pinned}/libxc-${Libxc_pinned}.tar.gz
      https://gitlab.com/libxc/libxc/-/archive/release-6.2.0/libxc-release-6.2.0.tar.gz
    )

  enable_language(C)

  # compile 3rd derivative code
  set(DISABLE_KXC OFF CACHE BOOL "" FORCE)
  # compile 4th derivative code
  set(DISABLE_LXC OFF CACHE BOOL "" FORCE)
  # disable compilation of testing infrastructure
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  # only build shared libraries
  set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
  # enable arch-dependent flags if they are used for VeloxChem
  set(ENABLE_XHOST ${ENABLE_ARCH_FLAGS} CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(Libxc)

  # Provide an alias, so linking to Libxc looks the same regardless if it was
  # found on the system or if it was fetched at configuration
  add_library(Libxc::xc ALIAS xc)
endif()
