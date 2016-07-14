include_directories (
    ${CMAKE_BINARY_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${BLAS_INCLUDE_DIRS}
    )

link_directories (
    ${CMAKE_BINARY_DIR}/lib
    ${Boost_LIBRARY_DIRS}
    )
