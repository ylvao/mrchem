
include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mwreprlib)
include_directories (${CMAKE_SOURCE_DIR}/src/mrgrid)

set (MRCPP_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrcpp.x)
set (MRCHEM_DRIVER_DEBUG ON)

add_executable(mrcpp.x mrcpp.cpp)

target_link_libraries(mrcpp.x mrgrid mwrepr)
target_link_libraries(mrcpp.x ${Boost_LIBRARIES})
target_link_libraries(mrcpp.x ${BLAS_LIBRARIES})

