include(DefaultPathsMRCPP)
include_directories (${Top_SOURCE_DIR}/mrchem)
include_directories (${Top_SOURCE_DIR}/mrcpp)

find_package(Boost REQUIRED program_options)
find_package(Boost REQUIRED serialization)

set (DEFAULT_LIBS mrcpp mrchem 
	${Boost_LIBRARIES}
	${BLAS_LIBRARIES}
	)

