#include "catch.hpp"

#include <cmath>
#include <iostream>

#include "TelePrompter.h"

TEST_CASE("Greeting the world of testing", "[hello_test]")
{
    println(0, "Hello world");
    REQUIRE( true );
}

