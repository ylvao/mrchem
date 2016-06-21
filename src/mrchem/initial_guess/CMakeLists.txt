include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwanalytic)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwcore)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwtrees)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwoperators)
include_directories (${CMAKE_SOURCE_DIR}/src/mrcpp/mwbuilders)

add_library(initial_guess STATIC 
    AOBasis.cpp
    AOContraction.cpp
    Intgrl.cpp
    OrbitalExp.cpp
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lsdalton.py.in lsdalton.py)

install(FILES 
        ${CMAKE_CURRENT_BINARY_DIR}/lsdalton.py
        periodictable.py
        DESTINATION ${PYTHON_SITE_INSTALL_DIR}
	)
