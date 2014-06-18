include(DefaultPathsMRCPP)
set(GTEST_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external/google-test/include)
set(GTEST_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/external/google-test)
set(GTEST_LIBRARIES gtest gtest_main)

include_directories (
    ${PROJECT_SOURCE_DIR}/src/mwreprlib
    ${PROJECT_SOURCE_DIR}/src/mrgrid
    ${GTEST_INCLUDE_DIRS}
    )
link_directories(${GTEST_LIBRARY_DIRS})

set(TESTING_LIBS mwrepr mrgrid
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
	NodeIndexTest
	BoundingBoxTest
	TreeTest 
	GridTest 
	)

foreach(test ${tests})
	add_executable(${test} ${test}.cpp)
	target_link_libraries(${test} ${TESTING_LIBS})
    	add_dependencies(${test} external-modules)
endforeach()

