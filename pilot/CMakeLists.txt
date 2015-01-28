
include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mwreprlib)
include_directories (${CMAKE_SOURCE_DIR}/src/mrgrid)
include_directories (${CMAKE_SOURCE_DIR}/src/mwutils)

set (MRGRID_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrgrid.bin)
set (MRCHEM_DRIVER_DEBUG ON)

add_executable(mrgrid.bin mrgrid.cpp)

target_link_libraries(mrgrid.bin mwutils mrgrid mwrepr)
target_link_libraries(mrgrid.bin ${Boost_LIBRARIES})
target_link_libraries(mrgrid.bin ${BLAS_LIBRARIES})

