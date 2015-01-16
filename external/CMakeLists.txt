include(ExternalProject)
add_custom_target(external-modules)

add_custom_target(git-submodule-init
    COMMAND git submodule init
    COMMAND git submodule update
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

ExternalProject_Add(submodule-libgetkw
    DOWNLOAD_COMMAND git submodule update
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libgetkw
    BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/libgetkw
    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    -DCMAKE_BUILD_TYPE=Release
    )
add_dependencies(submodule-libgetkw git-submodule-init)
add_dependencies(external-modules submodule-libgetkw)

if (ENABLE_XCFUN)
    enable_language(Fortran)
    ExternalProject_Add(submodule-xcfun
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/xcfun
        BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/xcfun
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
        -DCMAKE_BUILD_TYPE=Release
        )
    add_dependencies(submodule-xcfun git-submodule-init)
    add_dependencies(external-modules submodule-xcfun)
endif()

# vim:et:sw=4:ts=4:
