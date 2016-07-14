include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(pilot)
add_subdirectory(tools EXCLUDE_FROM_ALL)

configure_file (
    "${PROJECT_SOURCE_DIR}/config.h.in"
    "${PROJECT_BINARY_DIR}/config.h"
    )
