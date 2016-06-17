include(ExternalProject)
add_custom_target(external-modules)

add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )

ExternalProject_Add(submodule-libgetkw
    DOWNLOAD_COMMAND git submodule update
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libgetkw
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=Release
    )
add_dependencies(submodule-libgetkw git-submodule-init)
add_dependencies(external-modules submodule-libgetkw)

option(ENABLE_XCFUN "Use XCFun library" ON)
if (ENABLE_XCFUN)
    enable_language(Fortran)
    ExternalProject_Add(submodule-xcfun
        DOWNLOAD_COMMAND git submodule update
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/xcfun
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_BUILD_TYPE=Release
        )
    add_dependencies(submodule-xcfun git-submodule-init)
    add_dependencies(external-modules submodule-xcfun)
endif()

include_directories(${PROJECT_SOURCE_DIR}/external/catch)

# vim:et:sw=4:ts=4:
