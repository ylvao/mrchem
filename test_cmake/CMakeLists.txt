cmake_minimum_required(VERSION 3.14)
cmake_policy(SET CMP0144 NEW)

project(TestLibXC)

add_library(Libxc::libxc STATIC IMPORTED)

# Set properties for the imported library
set_target_properties(Libxc::libxc PROPERTIES
  IMPORTED_LOCATION "/home/ylva/work/libxc-mrchem/libxc/libxc_install/lib/libxc.a" # Path to the library
  INTERFACE_INCLUDE_DIRECTORIES "/home/ylva/work/libxc-mrchem/libxc/libxc_install/include" # Path to include files
)


# Set the path to libxc
set(CMAKE_PREFIX_PATH "/home/ylva/work/libxc-mrchem/libxc/libxc_install")

# Try to find libXC
find_package(Libxc)

if(TARGET Libxc::libxc)
  message(STATUS "LIBXC! Found libXC")
else()
  message(STATUS "LIBXC! Suitable libXC could not be located.")
endif()
