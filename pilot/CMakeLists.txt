project(Pilot)

include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src/mwreprlib)
include_directories (${CMAKE_SOURCE_DIR}/src/mrgrid)

find_package(PythonInterp REQUIRED)

set (MRGRID_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrgrid.bin)
set (MRCHEM_DRIVER_DEBUG ON)

add_executable(mrgrid.bin mrgrid.cpp)
add_dependencies(mrgrid.bin external-modules)
target_link_libraries(mrgrid.bin mrgrid mwrepr xcfun getkw)
target_link_libraries(mrgrid.bin ${Boost_LIBRARIES})
target_link_libraries(mrgrid.bin ${BLAS_LIBRARIES})

configure_file(${CMAKE_SOURCE_DIR}/src/mrgrid/mrgrid.in
    ${CMAKE_BINARY_DIR}/pilot/mrgrid)
execute_process(COMMAND
    chmod 755 ${CMAKE_BINARY_DIR}/pilot/mrgrid OUTPUT_QUIET)

