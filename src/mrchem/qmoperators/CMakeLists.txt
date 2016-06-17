include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(qmoperators STATIC 
    FockOperator.cpp
    KineticOperator.cpp
    MomentumOperator.cpp
    NuclearPotential.cpp
    Potential.cpp
    QMOperator.cpp
    XCFunctional.cpp
)

target_link_libraries(qmoperators
    mrchem
    chemistry
    initial_guess
    qmfunctions
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
    xcfun
)

install(TARGETS qmoperators
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

