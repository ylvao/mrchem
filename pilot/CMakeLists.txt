include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwarithmetics)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)

    target_link_libraries(mrchem-pilot.x mrchem mrcpp)
    target_link_libraries(mrchem-pilot.x ${Boost_LIBRARIES})
    target_link_libraries(mrchem-pilot.x ${BLAS_LIBRARIES})
endif()

