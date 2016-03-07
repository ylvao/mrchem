include_directories (${CMAKE_INSTALL_PREFIX}/include)

set (mrchem_SOURCES
        MREnv.cpp
        MREnv.h
	)

add_library(mrchem STATIC 
        ${mrchem_SOURCES}
        )

target_link_libraries(mrchem mrcpp xcfun)
add_dependencies(mrchem external-modules)

install(TARGETS mrchem
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
