include(${CMAKE_CURRENT_LIST_DIR}/testing_macros.cmake)

enable_testing()
include(CTest)
add_subdirectory(tests) # This must come last!!
