
set(mwrepr_SOURCES 
	constants.h
	eigen_disable_warnings.h
	macros.h
	mwrepr_declarations.h
	MathUtils.cpp
	MathUtils.h
	NodeBox.cpp
	NodeBox.h
	BoundingBox.cpp
	BoundingBox.h
	FunctionTree.cpp
	FunctionTree.h
	MRTree.cpp
	MRTree.h
	MWTree.cpp
	MWTree.h
	MRNode.cpp
	MRNode.h
	MWNode.cpp
	MWNode.h
	MWFilter.cpp
	MWFilter.h
	TreeProjector.h
	AnalyticProjector.cpp
	AnalyticProjector.h
	TreeAdaptor.h
	MWAdaptor.cpp
	MWAdaptor.h
	TreeBuilder.cpp
	TreeBuilder.h
	FunctionProjector.h
	FilterCache.cpp
	FilterCache.h
	FunctionNode.cpp
	FunctionNode.h
	ProjectedNode.cpp
	ProjectedNode.h
	GenNode.cpp
	GenNode.h
	NodeIndex.h
	parallel.cpp
	parallel.h
	TelePrompter.cpp
	TelePrompter.h
	TreeIterator.h
	HilbertPath.cpp
	HilbertPath.h
	HilbertIterator.h
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
        ScalingCache.h
        ScalingBasis.h
        LegendreBasis.h
        InterpolatingBasis.h
	)

add_library(mwrepr STATIC ${mwrepr_SOURCES})

install(TARGETS mwrepr 
	ARCHIVE DESTINATION lib
	CONFIGURATIONS Debug-devel Release-devel
	)

