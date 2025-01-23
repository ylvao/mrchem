#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MRCPP::mrcpp" for configuration "Release"
set_property(TARGET MRCPP::mrcpp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(MRCPP::mrcpp PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmrcpp.so.1"
  IMPORTED_SONAME_RELEASE "libmrcpp.so.1"
  )

list(APPEND _cmake_import_check_targets MRCPP::mrcpp )
list(APPEND _cmake_import_check_files_for_MRCPP::mrcpp "${_IMPORT_PREFIX}/lib/libmrcpp.so.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
