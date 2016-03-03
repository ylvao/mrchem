include(DefaultPathsMRCPP)
include_directories (${Top_SOURCE_DIR}/mwreprlib)
include_directories (${Top_SOURCE_DIR}/mwoperlib)

find_package(Boost REQUIRED program_options)
find_package(Boost REQUIRED serialization)

set (DEFAULT_LIBS mwrepr mwoper 
	${Boost_LIBRARIES}
	${BLAS_LIBRARIES}
	)

add_executable(ls_tree ls_tree.cpp)
target_link_libraries(ls_tree ${DEFAULT_LIBS})

add_executable(diff_trees diff_trees.cpp)
target_link_libraries(diff_trees ${DEFAULT_LIBS})

