include_directories (${CMAKE_INSTALL_PREFIX}/include)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)
link_directories (${Boost_LIBRARY_DIRS})

set (mrgridlib_SOURCES
	Atom.h
	AtomicElement.h
	Molecule.h
	PeriodicTable.cpp
	PeriodicTable.h
	MolecularGridGenerator.cpp
	MolecularGridGenerator.h
        MREnv.cpp
        MREnv.h
        mrgrid.in
	)

add_library(mrgrid STATIC 
        ${mrgridlib_SOURCES}
        )

target_link_libraries(mrgrid mwrepr)
add_dependencies(mrgrid external-modules)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mrgrid.in mrgrid)

install(TARGETS mrgrid
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
