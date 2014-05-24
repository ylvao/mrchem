add_library(mrgrid STATIC 
	Atom.h
	AtomicElement.h
	Molecule.h
	PeriodicTable.cpp
	PeriodicTable.h
	MolecularGridGenerator.cpp
	MolecularGridGenerator.h
	)

install(TARGETS mrgrid 
	ARCHIVE DESTINATION lib
	CONFIGURATIONS Debug-devel Release-devel
	)

