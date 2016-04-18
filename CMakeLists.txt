cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(MRCHEM CXX C)
enable_language(CXX)

set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 0)
set (PROJECT_VERSION_PATCH 0)
set (PROJECT_VERSION 
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})

option(ENABLE_XCFUN "Use XCFun library" OFF)
option(ENABLE_CUDA "Enable CUDA GPU acceleration support" OFF)
option(ENABLE_CODE_COVERAGE "Compile with code coverage turned on" OFF)

include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigCompilerFlags)

find_package(Eigen3 REQUIRED)
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost REQUIRED serialization system timer)

option(ENABLE_MPI "Enable MPI parallelization" OFF)
if (ENABLE_MPI)
    find_package(MPI)
    if (MPI_FOUND)
        find_package(Boost REQUIRED mpi)
        set(HAVE_MPI TRUE)
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif ()

option(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(HAVE_OPENMP TRUE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif ()

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

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
add_subdirectory(tools EXCLUDE_FROM_ALL)

option(MRCPP_INSTALL_FILTERS "Install MW filters" OFF)
if (MRCPP_INSTALL_FILTERS)
    install(DIRECTORY ${MW_FILTER_DIR}
        DESTINATION share
        )
    set (MW_FILTER_DIR ${CMAKE_INSTALL_PREFIX}/share/MultiwaveletFilters)
endif ()

option(ENABLE_TESTS "Enable test suite" ON)
if (ENABLE_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(tests)
endif()

find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_subdirectory(doc EXCLUDE_FROM_ALL)
endif()

configure_file (
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )

find_package(Git REQUIRED)

# Configure CPack
include(ConfigPackaging) 

include(FeatureSummary)
feature_summary(
   DESCRIPTION "Found packages:" 
   WHAT PACKAGES_FOUND
   )
