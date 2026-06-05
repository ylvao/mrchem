find_package(XCFun 2.1 CONFIG QUIET
  NO_CMAKE_PATH
  NO_CMAKE_PACKAGE_REGISTRY
  NO_CMAKE_SYSTEM_PACKAGE_REGISTRY
  )
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
else()
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)
  FetchContent_Declare(xcfun_sources
    QUIET
    URL
      https://github.com/dftlibs/xcfun/archive/v2.1.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

  set(CMAKE_BUILD_TYPE Release)
  set(ENABLE_TESTALL FALSE CACHE BOOL "" FORCE)
  set(XCFUN_MAX_ORDER 3)
  set(XCFUN_PYTHON_INTERFACE FALSE CACHE BOOL "" FORCE)

  # Remove this line to restore the "old" pbe behaviour when using XCFun
  add_compile_definitions(XCFUN_REF_PBEX_MU)

  FetchContent_MakeAvailable(xcfun_sources)

  if(TARGET xcfun)
    # Suppress all warnings from XCFun's own compilation (-w overrides any -Wall/-Wextra from CMAKE_CXX_FLAGS)
    target_compile_options(xcfun PRIVATE -w)
  endif()
endif()
