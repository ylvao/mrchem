include_directories (${CMAKE_INSTALL_PREFIX}/include)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)
link_directories (${Boost_LIBRARY_DIRS})

add_executable(mrchem.x mrchem.cpp)
add_dependencies(mrchem.x external-modules)

add_library(mrchem STATIC 
        mrchem.cpp
        MREnv.cpp
        MREnv.h
        )

target_link_libraries(mrchem mrcpp getkw)
add_dependencies(mrchem external-modules)

target_link_libraries(mrchem.x mrchem mrcpp getkw
	${Boost_LIBRARIES} ${BLAS_LIBRARIES}
        )


option (MRCHEM_DRIVER_DEBUG "Debug the driver (extmod)" OFF)

set (MRCHEM_EXECUTABLE ${CMAKE_INSTALL_PREFIX}/bin/mrchem.x)

install(TARGETS mrchem.x
	RUNTIME DESTINATION bin
	)
