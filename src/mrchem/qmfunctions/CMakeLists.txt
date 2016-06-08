include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(qmfunctions STATIC 
    Orbital.cpp
    OrbitalVector.cpp
)

target_link_libraries(qmfunctions
    mrchem
    chemistry
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
)

install(TARGETS qmfunctions
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

