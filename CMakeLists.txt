cmake_minimum_required(VERSION 2.8)

project(MRCHEM)

enable_language(C)
enable_language(CXX)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/custom
    )

include(code_coverage)
include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigCompilerFlags)

find_package(PythonInterp REQUIRED)
set(PYTHON_SITE_INSTALL_DIR lib/python2.7/site-packages)

find_package(Eigen3 REQUIRED)
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost REQUIRED serialization system timer)

include(mpi)
include(omp)
include(blas)

option(ENABLE_BLAS "Use BLAS backend for linear algebra" ON)
if (ENABLE_BLAS)
    if (BLAS_TYPE)
        find_package(BLAS COMPONENTS ${BLAS_TYPE})
    else()
        find_package(BLAS)
    endif()
endif ()

include_directories(${PROJECT_BINARY_DIR})

if (NOT DEFINED ${MW_FILTER_DIR})
    set (MW_FILTER_DIR "/home/jonas/dev/MultiwaveletFilters"
        CACHE STRING "Path to MW filters and cross-correlation coefs")
endif()

include(ConfigGitRevision)

if (NOT INITIAL_GUESS_PROGRAM)
    set (INITIAL_GUESS_PROGRAM "lsdalton")
endif()

if (NOT INITIAL_GUESS_PROGRAM_DIR)
    set (INITIAL_GUESS_PROGRAM_DIR "/unknown")
endif()

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
add_subdirectory(tools EXCLUDE_FROM_ALL)

include(mw_filters)
include(doxygen)
include(tests)

configure_file (
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )

include(feature_summary)
