include_directories (
    ${CMAKE_BINARY_DIR}/external/include
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${BLAS_INCLUDE_DIRS}
    )

link_directories (
    ${CMAKE_BINARY_DIR}/external/lib
    ${Boost_LIBRARY_DIRS}
    )
