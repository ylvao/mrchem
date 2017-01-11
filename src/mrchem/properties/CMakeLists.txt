include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(properties STATIC 
)

target_link_libraries(properties
    mrchem
    chemistry
    initial_guess
    qmfunctions
    qmoperators
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
    mwbuilders
    xcfun
)

install(TARGETS properties
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

