cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

project(MRCPP CXX C)
enable_language(CXX)

set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 0)
set (PROJECT_VERSION_PATCH 0)
set (PROJECT_VERSION 
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})


option(MRCPP_INSTALL_FILTERS "Install MW filters" OFF)
option(ENABLE_TESTS "Enable test suite" OFF)
option(ENABLE_XCFUN "Use XCFun library" ON)
option(ENABLE_BLAS "Use BLAS backend for linear algebra" ON)
option(ENABLE_MPI "Enable MPI parallelization" OFF)
option(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ENABLE_CUDA "Enable CUDA GPU acceleration support" OFF)
option(ENABLE_CODE_COVERAGE "Compile with code coverage turned on" OFF)

include(GNUInstallDirs)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigCompilerFlags)

#find_package(Python REQUIRED)
#set(PYTHON_SITE_INSTALL_DIR
#        lib/python${PYTHON_VERSION}/site-packages)

find_package(Eigen3 REQUIRED)
set(Boost_USE_MULTITHREADED OFF)
find_package(Boost REQUIRED serialization)

if (ENABLE_MPI)
    find_package(MPI)
    if (MPI_FOUND)
        find_package(Boost REQUIRED mpi)
        set(HAVE_MPI TRUE)
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif ()

if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(HAVE_OPENMP TRUE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif ()

if (ENABLE_CUDA)
    find_package(CUDA)
endif ()

if (ENABLE_BLAS)
    if (BLAS_TYPE)
        find_package(BLAS COMPONENTS ${BLAS_TYPE})
    else()
        find_package(BLAS)
    endif()
endif ()

include_directories(${PROJECT_BINARY_DIR})

if (NOT DEFINED ${MW_FILTER_DIR})
    set (MW_FILTER_DIR C"/home/jonas/dev/MultiwaveletFilters"
        CACHE STRING "Path to MW filters and cross-correlation coefs")
endif()

include(ConfigGitRevision)

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
add_subdirectory(tools EXCLUDE_FROM_ALL)

if (MRCPP_INSTALL_FILTERS)
    install(DIRECTORY ${MW_FILTER_DIR}
        DESTINATION share
        )
    set (MW_FILTER_DIR ${CMAKE_INSTALL_PREFIX}/share/MultiwaveletFilters)
endif ()

if (ENABLE_TESTS)
    set (PROJECT_REPOSITORY "git@repo.ctcc.no:mrcpp.git")
    include(ConfigTesting)
    add_subdirectory(test)
endif()

# Copy the examples in place for convenience 
if (NOT EXISTS ${CMAKE_BINARY_DIR}/mrgrid.inp)
    #set(_dir ${CMAKE_SOURCE_DIR}/examples)
    #file(COPY ${_dir}/mrchem.inp ${_dir}
        #DESTINATION ${CMAKE_BINARY_DIR}
        #)
    #unset(_dir)
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
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "MRCPP: Multi-Resoluton Chemistry Program Package")
set(CPACK_PACKAGE_VENDOR "CTCC")
set(CPACK_PACKAGE_CONTACT "Jonas Juselius <jonas.juselius@uit.no>")
set(CPACK_PACKAGE_EXECUTABLES "mrcpp" "The MRCPP program")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/INSTALL.rst")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_GENERATOR TGZ RPM DEB)
include(ConfigPackaging) 

include(FeatureSummary)
feature_summary(
   DESCRIPTION "Found packages:" 
   WHAT PACKAGES_FOUND
   )

