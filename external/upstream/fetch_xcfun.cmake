find_package(XCFun CONFIG QUIET)
if(TARGET XCFun::xcfun)
  get_property(_loc TARGET XCFun::xcfun PROPERTY LOCATION)
  message(STATUS "Found XCFun: ${_loc} (found version ${XCFun_VERSION})")
else()
  # FIXME XCFun needs the C compiler for the moment being
  #       Remove when updating to a latest and greatest that removes this annoyance
  enable_language(C)
  message(STATUS "Suitable XCFun could not be located. Fetching and building!")
  include(FetchContent)

  FetchContent_Populate(xcfun_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/dftlibs/xcfun.git
    GIT_TAG
      a486a3f1483925c84590cf7ea012cf177f67ae9b # Preferable to have a tag for a release
    GIT_SHALLOW
      1
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DENABLE_FC_SUPPORT=FALSE
      -DENABLE_TESTALL=TRUE
    CMAKE_CACHE_ARGS
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    )

  add_subdirectory(
    ${xcfun_sources_SOURCE_DIR}
    ${xcfun_sources_BINARY_DIR}
    )
endif()
