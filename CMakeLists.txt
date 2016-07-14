cmake_minimum_required(VERSION 2.8)

project(MRCHEM)

enable_language(C)
enable_language(CXX)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/custom
    )

include(GNUInstallDirs)



include(code_coverage)
include(ConfigProjectVersion)
include(ConfigSafeGuards)
include(ConfigCompilerFlags)

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



include(python)
include(eigen)
include(boost)
include(git)
include(initial_guess)
include(mw_filters)
include(doxygen)
include(tests)
include(main)
include(feature_summary)
