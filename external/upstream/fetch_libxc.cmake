find_package(Libxc 7.0 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET Libxc::xc)
  get_target_property(_loc Libxc::xc LOCATION)
  message(STATUS "Found Libxc: ${_loc} (found version ${Libxc_VERSION})")
else()
  message(STATUS "Suitable LibXC could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(libxc_sources
    QUIET
    URL
      https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "" FORCE)
  set(BUILD_TESTING     OFF CACHE BOOL "Build LibXC tests"          FORCE)
  set(ENABLE_TESTS      OFF CACHE BOOL "Enable LibXC tests"         FORCE)
  set(BUILD_SHARED_LIBS ON  CACHE BOOL "Build LibXC shared libs"    FORCE)
  set(ENABLE_FORTRAN    OFF CACHE BOOL "Build LibXC Fortran bindings" FORCE)
  set(DISABLE_FXC       ON  CACHE BOOL "Disable 2nd derivatives (Fxc)" FORCE)
  set(DISABLE_KXC       ON  CACHE BOOL "Disable 3rd derivatives (Kxc)" FORCE)
  set(DISABLE_LXC       ON  CACHE BOOL "Disable 4th derivatives (Lxc)" FORCE)

  FetchContent_MakeAvailable(libxc_sources)

  if(TARGET xc)
    if(NOT TARGET Libxc::xc)
      add_library(Libxc::xc ALIAS xc)
    endif()

    target_include_directories(xc
      INTERFACE
        $<BUILD_INTERFACE:${libxc_sources_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${libxc_sources_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/Libxc>
      )
  endif()
endif()
