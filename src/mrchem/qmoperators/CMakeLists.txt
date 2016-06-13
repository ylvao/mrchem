include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(qmoperators STATIC 
    KineticOperator.cpp
    MomentumOperator.cpp
    Potential.cpp
    FockOperator.cpp
    QMOperator.cpp
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
)

install(TARGETS qmoperators
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

