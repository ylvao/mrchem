include(DefaultPathsMRCPP)
include_directories (${Top_SOURCE_DIR}/mrchem)
include_directories (${Top_SOURCE_DIR}/mrcpp)


set (DEFAULT_LIBS mrcpp mrchem 
	${Boost_LIBRARIES}
	${BLAS_LIBRARIES}
	)

