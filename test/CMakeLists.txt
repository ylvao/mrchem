include(DefaultPathsMRCPP)
set(GTEST_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external/google-test/include)
set(GTEST_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/external/google-test)
set(GTEST_LIBRARIES gtest gtest_main)

include_directories (
    ${PROJECT_SOURCE_DIR}/src/mwreprlib
    ${PROJECT_SOURCE_DIR}/src/mwoperlib
    ${PROJECT_SOURCE_DIR}/src/mrchem
    ${GTEST_INCLUDE_DIRS}
    )
link_directories(${GTEST_LIBRARY_DIRS})

set(TESTING_LIBS mrchem mwrepr mwoper xcfun getkw
	${Boost_LIBRARIES}
	${BLAS_LIBRARIES}
    ${GTEST_LIBRARIES}
	)

if (MPI_FOUND)
    if (DEFINED ENV{MPIEXEC_PREFLAGS})
        separate_arguments(args UNIX_COMMAND $ENV{MPIEXEC_PREFLAGS})
        set (MPIEXEC_PREFLAGS ${args})
        unset(args)
    endif()
    if (DEFINED MPI_NUMPROC OR DEFINED $ENV{MPI_NUMPROC})    
        set (MPIRUN ${MPIEXEC} ${MPIEXEC_PREFLAGS} 
            ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROC} 
            )
    else()
        set (MPIRUN ${MPIEXEC} ${MPIEXEC_PREFLAGS})
    endif()
endif()

set(tests 
	TreeTest 
	SerializationTest 
	PoissonTest 
	TreeArithmeticTest 
	OperatorTest 
	XCFunTest 
	MoleculeTest
	SCFTest
	)

foreach(test ${tests})
	add_executable(${test} ${test}.cpp)
	target_link_libraries(${test} ${TESTING_LIBS})
    add_test(${test} ${MPIRUN} ${test} @mrchem.inp)
    add_dependencies(${test} external-modules)
endforeach()

add_custom_command(
	OUTPUT
		${CMAKE_CURRENT_BINARY_DIR}/@mrchem.inp 
		${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp 
		${CMAKE_CURRENT_BINARY_DIR}/test.bas
		${CMAKE_CURRENT_BINARY_DIR}/test.dens
	DEPENDS	
		${CMAKE_CURRENT_SOURCE_DIR}/mrchem.inp 
		${CMAKE_CURRENT_SOURCE_DIR}/test.bas
		${CMAKE_CURRENT_SOURCE_DIR}/test.dens
	WORKING_DIRECTORY 
		${CMAKE_CURRENT_BINARY_DIR}
	COMMAND 
		${CMAKE_COMMAND} -E copy 
			${CMAKE_CURRENT_SOURCE_DIR}/mrchem.inp ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND 
		${CMAKE_COMMAND} -E copy 
			${CMAKE_CURRENT_SOURCE_DIR}/test.bas ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND 
		${CMAKE_COMMAND} -E copy 
			${CMAKE_CURRENT_SOURCE_DIR}/test.dens ${CMAKE_CURRENT_BINARY_DIR}
	COMMAND 
		${CMAKE_BINARY_DIR}/mrchem -D -p 0
	)

add_custom_target(setup_test_data ALL
	DEPENDS	
	${CMAKE_CURRENT_BINARY_DIR}/mrchem.inp 
	${CMAKE_CURRENT_BINARY_DIR}/test.bas
	${CMAKE_CURRENT_BINARY_DIR}/test.dens
	)
add_dependencies(setup_test_data submodule-libgetkw)
