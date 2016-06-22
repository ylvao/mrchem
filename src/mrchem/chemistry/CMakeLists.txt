include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(chemistry STATIC 
    Molecule.cpp
    PeriodicTable.cpp
)

target_link_libraries(chemistry
    mrchem
    qmfunctions
    properties
    scf_solver
    mrcpp
    mwcore
    mwtrees
    mwoperators
    mwbuilders
)

install(TARGETS chemistry
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

