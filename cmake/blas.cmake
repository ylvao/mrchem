option(ENABLE_BLAS "Use BLAS backend for linear algebra" ON)

if(ENABLE_BLAS)
    if(BLAS_TYPE)
        find_package(BLAS COMPONENTS ${BLAS_TYPE})
    else()
        find_package(BLAS)
    endif()
endif()
