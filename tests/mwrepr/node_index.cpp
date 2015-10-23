#include "catch.hpp"

#include "NodeIndex.h"

// Factory helper functions
template<int D> void initialize(NodeIndex<D> **idx);
template<int D> void finalize(NodeIndex<D> **idx);
template<int D> void testInitial(const NodeIndex<D> *idx);

TEST_CASE("NodeIndex constructor 1D", "[node_index_constructor_1D], [node_index]") {
    NodeIndex<1> *nIdx = 0;
    initialize<1>(&nIdx);

    SECTION("Constructor") {
        testInitial<1>(nIdx);
    }

    SECTION("Copy constructor") {
        NodeIndex<1> *cIdx = new NodeIndex<1>(*nIdx);
        testInitial<1>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = 1;
        NodeIndex<1> *cIdx = new NodeIndex<1>(*nIdx, i);
        REQUIRE( cIdx->getScale() == (nIdx->getScale() + 1) );
        REQUIRE( cIdx->getRankId() == nIdx->getRankId() );
        finalize(&cIdx);
    }

    SECTION("Default constructor") {
        NodeIndex<1> *cIdx = new NodeIndex<1>();
        REQUIRE( cIdx->getRankId() < 0 );
        SECTION("Assignment operator") {
            *cIdx = *nIdx;
            testInitial<1>(cIdx);
        }
        finalize(&cIdx);
    }
    finalize(&nIdx);
}

TEST_CASE("NodeIndex constructor 2D", "[node_index_constructor_2D], [node_index]") {
    NodeIndex<2> *nIdx = 0;
    initialize<2>(&nIdx);

    SECTION("Constructor") {
        testInitial<2>(nIdx);
    }

    SECTION("Copy constructor") {
        NodeIndex<2> *cIdx = new NodeIndex<2>(*nIdx);
        testInitial<2>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = 2;
        NodeIndex<2> *cIdx = new NodeIndex<2>(*nIdx, i);
        REQUIRE( cIdx->getScale() == (nIdx->getScale() + 1) );
        REQUIRE( cIdx->getRankId() == nIdx->getRankId() );
        finalize(&cIdx);
    }

    SECTION("Default constructor") {
        NodeIndex<2> *cIdx = new NodeIndex<2>();
        REQUIRE( cIdx->getRankId() < 0 );
        SECTION("Assignment operator") {
            *cIdx = *nIdx;
            testInitial<2>(cIdx);
        }
        finalize(&cIdx);
    }
    finalize(&nIdx);
}

TEST_CASE("NodeIndex constructor 3D", "[node_index_constructor_3D], [node_index]") {
    NodeIndex<3> *nIdx = 0;
    initialize<3>(&nIdx);

    SECTION("Constructor") {
        testInitial<3>(nIdx);
    }

    SECTION("Copy constructor") {
        NodeIndex<3> *cIdx = new NodeIndex<3>(*nIdx);
        testInitial<3>(cIdx);
        finalize(&cIdx);
    }

    SECTION("Child constructor") {
        int i = 4;
        NodeIndex<3> *cIdx = new NodeIndex<3>(*nIdx, i);
        REQUIRE( cIdx->getScale() == (nIdx->getScale() + 1) );
        REQUIRE( cIdx->getRankId() == nIdx->getRankId() );
        finalize(&cIdx);
    }

    SECTION("Default constructor") {
        NodeIndex<3> *cIdx = new NodeIndex<3>();
        REQUIRE( cIdx->getRankId() < 0 );
        SECTION("Assignment operator") {
            *cIdx = *nIdx;
            testInitial<3>(cIdx);
        }
        finalize(&cIdx);
    }
    finalize(&nIdx);
}

SCENARIO("Node indices can be compared", "[node_index_compare], [node_index]") {
    GIVEN("Two identical node indices aIdx and bIdx in 1D") {
        NodeIndex<1> *aIdx = 0;
        NodeIndex<1> *bIdx = 0;
        initialize<1>(&aIdx);
        initialize<1>(&bIdx);
        THEN("aIdx == bIdx") {
            REQUIRE( *aIdx == *bIdx );
            REQUIRE_FALSE( *aIdx != *bIdx );
        }
        WHEN("aIdx is given a deeper scale") {
            aIdx->setScale(2);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("bIdx is given a coarser scale") {
            bIdx->setScale(-10);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("aIdx is given a different rank") {
            aIdx->setRankId(5);
            THEN("aIdx == bIdx") {
                REQUIRE( *aIdx == *bIdx );
                REQUIRE_FALSE( *aIdx != *bIdx );
            }
        }
        WHEN("aIdx is given a different translation") {
            int l[1] = {5};
            aIdx->setTranslation(l);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        finalize<1>(&aIdx);
        finalize<1>(&bIdx);
    }

    GIVEN("Two identical node indices aIdx and bIdx in 2D") {
        NodeIndex<2> *aIdx = 0;
        NodeIndex<2> *bIdx = 0;
        initialize<2>(&aIdx);
        initialize<2>(&bIdx);
        THEN("aIdx == bIdx") {
            REQUIRE( *aIdx == *bIdx );
            REQUIRE_FALSE( *aIdx != *bIdx );
        }
        WHEN("aIdx is given a deeper scale") {
            aIdx->setScale(2);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("bIdx is given a coarser scale") {
            bIdx->setScale(-10);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("aIdx is given a different rank") {
            aIdx->setRankId(5);
            THEN("aIdx == bIdx") {
                REQUIRE( *aIdx == *bIdx );
                REQUIRE_FALSE( *aIdx != *bIdx );
            }
        }
        WHEN("aIdx is given a different translation") {
            int l[2] = {5, 0};
            aIdx->setTranslation(l);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        finalize<2>(&aIdx);
        finalize<2>(&bIdx);
    }

    GIVEN("Two identical node indices aIdx and bIdx in 3D") {
        NodeIndex<3> *aIdx = 0;
        NodeIndex<3> *bIdx = 0;
        initialize<3>(&aIdx);
        initialize<3>(&bIdx);
        THEN("aIdx == bIdx") {
            REQUIRE( *aIdx == *bIdx );
            REQUIRE_FALSE( *aIdx != *bIdx );
        }
        WHEN("aIdx is given a deeper scale") {
            aIdx->setScale(2);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("bIdx is given a coarser scale") {
            bIdx->setScale(-10);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        WHEN("aIdx is given a different rank") {
            aIdx->setRankId(5);
            THEN("aIdx == bIdx") {
                REQUIRE( *aIdx == *bIdx );
                REQUIRE_FALSE( *aIdx != *bIdx );
            }
        }
        WHEN("aIdx is given a different translation") {
            int l[3] = {-1, 2, 1};
            aIdx->setTranslation(l);
            THEN("aIdx != bIdx") {
                REQUIRE( *aIdx != *bIdx );
                REQUIRE_FALSE( *aIdx == *bIdx );
            }
        }
        finalize<3>(&aIdx);
        finalize<3>(&bIdx);
    }
}

template<int D> void initialize(NodeIndex<D> **idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");
    if (*idx != 0) MSG_FATAL("Invalid argument");
    int scale = 1;
    int rank = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = d-1;
    }
    *idx = new NodeIndex<D>(scale, l, rank);
}

template<int D> void finalize(NodeIndex<D> **idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");
    if (*idx == 0) MSG_FATAL("Invalid argument");
    delete *idx;
    *idx = 0;
}

template<int D> void testInitial(const NodeIndex<D> *idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");

    const int scale = 1;
    REQUIRE( scale == idx->getScale() );

    const int rank = 1;
    REQUIRE( rank == idx->getRankId() );

    for (int d = 0; d < D; d++) {
        const int l = d-1;
        REQUIRE( l == idx->getTranslation(d) );
        REQUIRE( l == idx->getTranslation()[d] );
    }
}

