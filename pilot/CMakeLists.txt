include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem)

set (MRCHEM_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrchem.x)
set (MRCHEM_DRIVER_DEBUG ON)

add_executable(mrchem.x mrchem.cpp)

target_link_libraries(mrchem.x mrchem mrcpp)
target_link_libraries(mrchem.x ${Boost_LIBRARIES})
target_link_libraries(mrchem.x ${BLAS_LIBRARIES})

