set(mwutils_SOURCES 
    AOBasis.cpp
    AOBasis.h
    AOContraction.cpp
    AOContraction.h
    DensExp.cpp
    DensExp.h
    OrbExp.cpp
    OrbExp.h
    #Density.cpp
    #Density.h
    GaussExp.cpp
    GaussExp.h
    GaussFunc.cpp
    GaussFunc.h
    GaussPoly.cpp
    GaussPoly.h
    Gaussian.cpp
    Gaussian.h
    Intgrl.cpp
    Intgrl.h
	)

add_library(mwutils STATIC ${mwutils_SOURCES})

install(TARGETS mwutils 
	ARCHIVE DESTINATION lib
	CONFIGURATIONS Debug-devel Release-devel
	)


target_link_libraries(mwutils mrgrid mwrepr)
add_dependencies(mwutils external-modules)

install(TARGETS mwutils
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
