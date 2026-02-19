# include(FetchContent)

# # 1. Essential: Make sure Fortran is enabled at the top level
# # If this isn't in your main CMakeLists.txt, the fetch will fail to create targets.
# get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
# if(NOT "Fortran" IN_LIST languages)
#     message(FATAL_ERROR "LibDftD3 requires Fortran. Please add Fortran to your project() call in the main CMakeLists.txt.")
# endif()

# find_package(LibDftD3 CONFIG QUIET)

# if(NOT TARGET LibDftD3::libdftd3)
#   # Use lowercase consistently
#   FetchContent_Declare(
#     libdftd3
#     GIT_REPOSITORY https://github.com/cuanto/libdftd3.git
#     GIT_TAG        master
#   )
  
#   # Manual population to ensure we have the source directory variable
#   FetchContent_GetProperties(libdftd3)
#   if(NOT libdftd3_POPULATED)
#     FetchContent_Populate(libdftd3)
#     add_subdirectory(${libdftd3_SOURCE_DIR} ${libdftd3_BINARY_DIR})
#   endif()

#   # Based on the 'cuanto/libdftd3' repo, the target name is exactly 'dftd3'
#   if(TARGET dftd3)
#       add_library(LibDftD3::libdftd3 ALIAS dftd3)
#       message(STATUS "SUCCESS: Aliased dftd3 to LibDftD3::libdftd3")
#   else()
#       # If we reach here, the subdirectory was added but 'dftd3' doesn't exist.
#       message(FATAL_ERROR "Fetched LibDftD3 but target 'dftd3' was not found in ${libdftd3_SOURCE_DIR}")
#   endif()
# endif()