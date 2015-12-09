include_directories (${CMAKE_INSTALL_PREFIX}/include)

set (mrgridlib_SOURCES
	Atom.cpp
	Atom.h
	AtomicElement.h
	Molecule.cpp
	Molecule.h
	PeriodicTable.cpp
	PeriodicTable.h
	#MolecularGridGenerator.cpp
	#MolecularGridGenerator.h
        MREnv.cpp
        MREnv.h
        #XCFunctional.cpp
        #XCFunctional.h
        #XCProjector.cpp
        #XCProjector.h
	)

add_library(mrgrid STATIC 
        ${mrgridlib_SOURCES}
        )

target_link_libraries(mrgrid mwrepr xcfun)
add_dependencies(mrgrid external-modules)

install(TARGETS mrgrid
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
