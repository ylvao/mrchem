#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "XCFun::xcfun" for configuration "Release"
set_property(TARGET XCFun::xcfun APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(XCFun::xcfun PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libxcfun.so.2"
  IMPORTED_SONAME_RELEASE "libxcfun.so.2"
  )

list(APPEND _cmake_import_check_targets XCFun::xcfun )
list(APPEND _cmake_import_check_files_for_XCFun::xcfun "${_IMPORT_PREFIX}/lib/libxcfun.so.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
