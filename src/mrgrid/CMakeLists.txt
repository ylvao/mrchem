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
	)

if (XCFUN_FOUND)
        set (mrgridlib_sources ${mrgridlib_SOURCES} XCFun.cpp)
        set (mrgridlib_sources ${mrgridlib_SOURCES} XCFun.h)
endif()

add_library(mrgrid STATIC 
        ${mrgridlib_SOURCES}
        )

target_link_libraries(mrgrid mwrepr xcfun getkw)
add_dependencies(mrgrid external-modules)

install(TARGETS mrgrid
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
