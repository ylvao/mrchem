
set(mwrepr_SOURCES 
	checks.h
	constants.h
	eigen_disable_warnings.h
	macros.h
	mwrepr_declarations.h
	MathUtils.cpp
	MathUtils.h
	GridNode.cpp
	GridNode.h
	GridNodeBox.cpp
	GridNodeBox.h
	GridGenerator.cpp
	GridGenerator.h
	MRGrid.cpp
	MRGrid.h
	MultiException.cpp
	MultiException.h
	NodeIndex.h
	parallel.cpp
	parallel.h
	TelePrompter.cpp
	TelePrompter.h
	TreeIterator.h
	LebesgueIterator.h
	GaussQuadrature.cpp
	GaussQuadrature.h
	QuadratureCache.cpp
	QuadratureCache.h
	RepresentableFunction.cpp
	RepresentableFunction.h
	ObjectCache.h
	LegendrePoly.cpp
	LegendrePoly.h
	Polynomial.cpp
	Polynomial.h
	)

add_library(mwrepr STATIC ${mwrepr_SOURCES})

install(TARGETS mwrepr 
	ARCHIVE DESTINATION lib
	CONFIGURATIONS Debug-devel Release-devel
	)

