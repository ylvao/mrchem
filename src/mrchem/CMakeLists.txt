include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

include_directories (chemistry)

add_subdirectory (chemistry)

add_executable(mrchem.x mrchem.cpp)
add_dependencies(mrchem.x external-modules)

add_library(mrchem STATIC 
    mrchem.cpp
    MREnv.cpp
)

target_link_libraries(mrchem 
    chemistry 
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
)

target_link_libraries(mrchem.x 
    chemistry 
    mrchem 
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
)

add_dependencies(mrchem external-modules)

option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)

set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)

install(TARGETS mrchem.x
	RUNTIME DESTINATION bin
	)
