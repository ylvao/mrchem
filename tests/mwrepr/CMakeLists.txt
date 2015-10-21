set_property(GLOBAL APPEND PROPERTY TestSources ${CMAKE_CURRENT_LIST_DIR}/hello_test.cpp)
add_test(NAME hello_test COMMAND unit_tests.x [hello_test])

