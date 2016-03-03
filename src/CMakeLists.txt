project(Top CXX)

include(DefaultPathsMRCPP)
include_directories(mwreprlib)
include_directories(mwoperlib)
include_directories(mrchem)

add_subdirectory(mwreprlib)
add_subdirectory(mwoperlib)
add_subdirectory(mrchem)

