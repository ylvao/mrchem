
set(mwrepr_SOURCES 
	macros.h
	constants.h
	eigen_disable_warnings.h
	mwrepr_declarations.h
	parallel.cpp
	parallel.h
	MathUtils.cpp
	MathUtils.h
        MultiResolutionAnalysis.h
	NodeBox.cpp
	NodeBox.h
	BoundingBox.cpp
	BoundingBox.h
	MWFilter.cpp
	MWFilter.h
	FilterCache.cpp
	FilterCache.h
	HilbertPath.cpp
	HilbertPath.h
	NodeIndex.h
	TelePrompter.cpp
	TelePrompter.h
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
        MRTree.cpp
        MRTree.h
        MRNode.cpp
        MRNode.h
        MWTree.cpp
        MWTree.h
        MWNode.cpp
        MWNode.h
        FunctionTree.cpp
        FunctionTree.h
        FunctionNode.cpp
        FunctionNode.h
        ProjectedNode.cpp
        ProjectedNode.h
        GenNode.cpp
        GenNode.h
	)

add_library(mwrepr STATIC ${mwrepr_SOURCES})

install(TARGETS mwrepr 
	ARCHIVE DESTINATION lib
	CONFIGURATIONS Debug-devel Release-devel
	)

