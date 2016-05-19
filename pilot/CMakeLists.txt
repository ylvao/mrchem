include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mrchem)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrchem.cpp)
    add_executable(mrchem-pilot.x mrchem.cpp)

    target_link_libraries(mrchem-pilot.x 
        mrchem 
        mrcpp 
        mwanalytic 
        mwcore
        mwtrees 
        mwoperators 
        mwbuilders
        ${Boost_LIBRARIES}
        ${BLAS_LIBRARIES}
    )
endif()

