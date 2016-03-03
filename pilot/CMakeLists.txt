project(Pilot)

include(DefaultPathsMRCPP)
include_directories (${Top_SOURCE_DIR}/mwreprlib)
include_directories (${Top_SOURCE_DIR}/mwoperlib)
include_directories (${Top_SOURCE_DIR}/mrchem)

if (XCFUN_FOUND)
	set (mrchemlib_sources ${mrchemlib_sources} XCFun.cpp)
endif()

find_package(PythonInterp REQUIRED)

set (INST_LIBDIR ${CMAKE_BINARY_DIR}/external/share/libgetkw)
set (MRCHEM_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrchem-pilot.bin)
set (MRCHEM_DRIVER_DEBUG ON)

if (EXISTS ${Pilot_SOURCE_DIR}/mrchem.cpp)
	if (ENABLE_CUDA)
		cuda_add_executable(mrchem-pilot.bin mrchem.cpp CUDA_test.cu)
	else()
		add_executable(mrchem-pilot.bin mrchem.cpp)
	endif()
    add_dependencies(mrchem-pilot.bin external-modules)
	target_link_libraries(mrchem-pilot.bin mrchem mwoper mwrepr xcfun getkw)
	target_link_libraries(mrchem-pilot.bin ${Boost_LIBRARIES})
	target_link_libraries(mrchem-pilot.bin ${BLAS_LIBRARIES})
	configure_file(${CMAKE_SOURCE_DIR}/src/mrchem/mrchem.in 
		${CMAKE_BINARY_DIR}/pilot/mrchem)
	execute_process(COMMAND 
		chmod 755 ${CMAKE_BINARY_DIR}/pilot/mrchem OUTPUT_QUIET) 
endif()

set(pilot_inputs
	mrchem.inp 
	test.bas 
	test.mo
	test.en
	)

foreach (input_file ${pilot_inputs})
	if (NOT EXISTS ${Pilot_BINARY_DIR}/${input_file})
		file(COPY ${input_file} 
			DESTINATION ${Pilot_BINARY_DIR})
	endif()
endforeach()

