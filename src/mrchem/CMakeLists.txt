include_directories (${CMAKE_INSTALL_PREFIX}/include) 

include_directories (chemistry)
include_directories (initial_guess)

add_subdirectory (chemistry)
add_subdirectory (initial_guess)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)
link_directories (${Boost_LIBRARY_DIRS})

add_executable(mrchem.x mrchem.cpp)
add_dependencies(mrchem.x external-modules)

add_library(mrchem STATIC 
    mrchem.in
    mrchem.cpp
    MREnv.cpp
    SCFDriver.cpp
)

target_link_libraries(mrchem 
    chemistry 
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
    getkw
)

add_dependencies(mrchem external-modules)

target_link_libraries(mrchem.x 
    mrchem 
    chemistry
    mrcpp 
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
    getkw
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
)

option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)

set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrchem.in mrchem)

install(TARGETS mrchem.x
	RUNTIME DESTINATION bin
	)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/mrchem
	DESTINATION bin
	)
