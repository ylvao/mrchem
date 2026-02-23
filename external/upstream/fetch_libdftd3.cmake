include(FetchContent)

FetchContent_Declare(libdftd3
  GIT_REPOSITORY https://github.com/cuanto/libdftd3.git
  GIT_TAG        master
)

FetchContent_GetProperties(libdftd3)
if(NOT libdftd3_POPULATED)
  FetchContent_Populate(libdftd3)
endif()

# new libdftd3 version places Fortran sources under lib/
# verify that at least one expected file is present
if(NOT EXISTS "${libdftd3_SOURCE_DIR}/lib/common.f90")
    message(FATAL_ERROR "libdftd3 sources not available; check network/git")
endif()

# gather all Fortran source files from the lib directory
file(GLOB dftd3_sources
    "${libdftd3_SOURCE_DIR}/lib/*.f90"
)

# alternatively, could list them explicitly:
# set(dftd3_sources
#     "${libdftd3_SOURCE_DIR}/lib/common.f90"
#     "${libdftd3_SOURCE_DIR}/lib/sizes.f90"
#     "${libdftd3_SOURCE_DIR}/lib/pars.f90"
#     "${libdftd3_SOURCE_DIR}/lib/core.f90"
#     "${libdftd3_SOURCE_DIR}/lib/wrapper.f90"
#     "${libdftd3_SOURCE_DIR}/lib/api.f90"
#     # extras.f90 lives in prg/ not lib; may not be needed for library
# )

# 1. Create the library target (using a simple name)
  add_library(dftd3_lib STATIC ${dftd3_sources})
  
  # 2. IMPORTANT: Create the Alias that MRChem is looking for
  # This stops the linker from looking for a literal "-lLibDftD3::libdftd3"
  add_library(LibDftD3::libdftd3 ALIAS dftd3_lib)

  # 3. Set standard properties
  set_target_properties(dftd3_lib PROPERTIES 
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/fortran_modules"
  )
  
  target_include_directories(dftd3_lib PUBLIC "${libdftd3_SOURCE_DIR}/src")