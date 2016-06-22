include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(scf_solver STATIC 
    GroundStateSolver.cpp
    HelmholtzOperatorSet.cpp
    SCF.cpp
)

target_link_libraries(scf_solver
    mrchem
    chemistry
    initial_guess
    qmfunctions
    qmoperators
    properties
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
    xcfun
)

install(TARGETS scf_solver
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

