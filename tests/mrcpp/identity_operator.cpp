#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "IdentityOperator.h"
#include "IdentityKernel.h"
#include "CrossCorrelationGenerator.h"
#include "OperatorTree.h"
#include "BandWidth.h"

namespace identity_operator {

template<int D> void applyIdentity();

TEST_CASE("Initialize identity operator", "[init_identity], [identity_operator], [mw_operator]") {
    double exp_prec  = 1.0e-6;
    double proj_prec = 1.0e-6;
    double ccc_prec  = 1.0e-4;

    const int n = -6;
    const int k = 5;

    SECTION("Initialize identity kernel") {
        IdentityKernel id_kern(exp_prec);
        REQUIRE( id_kern.size() == 1 );

        SECTION("Project identity kernel") {
            int l = -1;
            int nbox = 2;
            NodeIndex<1> idx(n, &l);
            BoundingBox<1> box(idx, &nbox);

            InterpolatingBasis basis(2*k+1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);
            GridGenerator<1> G(kern_mra);
            MWProjector<1> Q(kern_mra, proj_prec);

            FunctionTree<1> *kern_tree = G(id_kern);
            Q(*kern_tree, id_kern);
            REQUIRE( kern_tree->integrate() == Approx(1.0).epsilon(proj_prec) );

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                CrossCorrelationGenerator G(oper_mra, ccc_prec);
                OperatorTree *oper_tree = G(*kern_tree);

                oper_tree->calcBandWidth(1.0);
                BandWidth bw_1 = oper_tree->getBandWidth();
                oper_tree->clearBandWidth();

                oper_tree->calcBandWidth(0.001);
                BandWidth bw_2 = oper_tree->getBandWidth();
                oper_tree->clearBandWidth();

                oper_tree->calcBandWidth(-1.0);
                BandWidth bw_3 = oper_tree->getBandWidth();
                oper_tree->clearBandWidth();

                for (int i = 0; i < oper_tree->getDepth(); i++) {
                    REQUIRE( bw_1.getMaxWidth(i) <= bw_2.getMaxWidth(i) );
                    REQUIRE( bw_2.getMaxWidth(i) <= bw_3.getMaxWidth(i) );
                }

                delete oper_tree;
            }
            delete kern_tree;
        }
    }
}

TEST_CASE("Apply identity operator", "[apply_identity], [identity_operator], [mw_operator]") {
    SECTION("1D") {
        applyIdentity<1>();
    }
    SECTION("2D") {
        applyIdentity<2>();
    }
    SECTION("3D") {
        applyIdentity<3>();
    }
}

template<int D> void applyIdentity() {
    double proj_prec = 1.0e-3;
    double apply_prec = 1.0e-3;
    double build_prec = 1.0e-4;

    MultiResolutionAnalysis<D> *mra = 0;
    GaussFunc<D> *fFunc = 0;
    initialize(&fFunc);
    initialize(&mra);

    MWProjector<D> Q(*mra, proj_prec);
    FunctionTree<D> *fTree = Q(*fFunc);

    GridGenerator<D> G(*mra);
    FunctionTree<D> *gTree = G();

    IdentityOperator<D> I(*mra, apply_prec, build_prec);
    I(*gTree, *fTree);

    REQUIRE( gTree->getDepth() <= fTree->getDepth() );
    REQUIRE( gTree->getNNodes() <= fTree->getNNodes() );
    REQUIRE( gTree->integrate() == Approx(fTree->integrate()).epsilon(apply_prec) );

    delete gTree;
    delete fTree;
    finalize(&fFunc);
    finalize(&mra);
}

} // namespace
