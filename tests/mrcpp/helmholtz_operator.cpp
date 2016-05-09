#include "catch.hpp"

#include "factory_functions.h"
#include "HelmholtzOperator.h"
#include "OperatorTreeVector.h"
#include "MWProjector.h"
#include "MWMultiplier.h"
#include "MWAdder.h"
#include "BandWidth.h"
#include "CrossCorrelationGenerator.h"
#include "HydrogenicFunction.h"
#include "NuclearFunction.h"

using namespace std;

namespace helmholtz_operator {

TEST_CASE("Helmholtz' kernel", "[init_helmholtz], [helmholtz_operator], [mw_operator]") {
    const double mu = 1.0;
    const double r_min = 1.0e-3;
    const double r_max = 1.0e+1;
    const double exp_prec  = 1.0e-4;
    const double proj_prec = 1.0e-3;
    const double ccc_prec  = 1.0e-3;
    const double band_prec  = 1.0e-3;

    const int n = -3;
    const int k = 5;

    SECTION("Initialize Helmholtz' kernel") {
        HelmholtzKernel helmholtz(mu, exp_prec, r_min, r_max);
        REQUIRE( helmholtz.size() == 24 );

        int foo = 0;
        double x = r_min;
        while (x < r_max) {
            REQUIRE( helmholtz.evalf(&x) == Approx(exp(-mu*x)/x).epsilon(exp_prec) );
            x *= 1.5;
        }
        SECTION("Project Helmholtz' kernel") {
            int l = -1;
            int nbox = 2;
            NodeIndex<1> idx(n, &l);
            BoundingBox<1> box(idx, &nbox);

            InterpolatingBasis basis(2*k+1);
            MultiResolutionAnalysis<1> kern_mra(box, basis);

            MWProjector<1> Q(kern_mra, proj_prec);

            FunctionTreeVector<1> kern_vec;
            for (int i = 0; i < helmholtz.size(); i++) {
                Gaussian<1> &kern_gauss = *helmholtz[i];
                FunctionTree<1> *kern_tree = Q(kern_gauss);
                kern_vec.push_back(*kern_tree);
            }

            SECTION("Build operator tree by cross correlation") {
                NodeIndex<2> idx(n);
                BoundingBox<2> box(idx);

                InterpolatingBasis basis(k);
                MultiResolutionAnalysis<2> oper_mra(box, basis);

                CrossCorrelationGenerator G(oper_mra, ccc_prec);

                OperatorTreeVector oper_vec;
                for (int i = 0; i < kern_vec.size(); i++) {
                    FunctionTree<1> &kern_tree = *kern_vec[i];
                    OperatorTree *oper_tree = G(kern_tree);
                    oper_vec.push_back(*oper_tree);

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
                }
                oper_vec.calcBandWidths(band_prec);
                REQUIRE( oper_vec.getMaxBandWidth(3) == 3 );
                REQUIRE( oper_vec.getMaxBandWidth(7) == 5 );
                REQUIRE( oper_vec.getMaxBandWidth(13) == 9 );
                REQUIRE( oper_vec.getMaxBandWidth(19) == -1 );

                for (int i = 0; i < oper_vec.size(); i++) {
                    delete oper_vec[i];
                }
                oper_vec.clear();
            }
            for (int i = 0; i < kern_vec.size(); i++) {
                delete kern_vec[i];
            }
            kern_vec.clear();
        }
    }
}

TEST_CASE("Apply Helmholtz' operator", "[apply_helmholtz], [helmholtz_operator], [mw_operator]") {
    double proj_prec = 1.0e-3;
    double apply_prec = 1.0e-2;
    double build_prec = 1.0e-3;

    // Computational domain [-32.0, 32.0]
    int scale = -5;
    int corner[3] = {-1, -1, -1};
    int nbox[3] = {2, 2, 2};
    NodeIndex<3> idx(scale, corner);
    BoundingBox<3> box(idx, nbox);

    int order = 5;
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(box, basis);

    FunctionTreeVector<3> tree_vec;
    MWAdder<3> add(MRA);
    MWMultiplier<3> mult(MRA);
    MWProjector<3> Q(MRA, proj_prec);
    GridGenerator<3> G(MRA);

    int n = 2;                  // Principal quantum number
    int l = 1;                  // Angular quantum number
    int m_l = 2;                // Magnetic quantum number
    double Z = 1.0;             // Nuclear charge
    double E = -Z/(2.0*n*n);    // Total energy

    double mu = sqrt(-2*E);
    HelmholtzOperator H(MRA, mu, apply_prec, build_prec);

    double pos[3] = {0.0, 0.0, 0.0};
    HydrogenicFunction hFunc(n, l, m_l, Z, pos);
    FunctionTree<3> *psi_n = Q(hFunc);

    NuclearFunction nucFunc;
    nucFunc.addNucleus(pos, Z, 1.0e-6);
    FunctionTree<3> *V = Q(nucFunc);

    tree_vec.push_back(*V);
    tree_vec.push_back(*psi_n);
    FunctionTree<3> *Vpsi = mult(tree_vec);
    tree_vec.clear();

    FunctionTree<3> *psi_np1 = G(*psi_n);
    H(*psi_np1, *Vpsi);
    *psi_np1 *= -1.0/(2.0*pi);

    double norm = sqrt(psi_np1->getSquareNorm());
    REQUIRE( norm == Approx(1.0).epsilon(apply_prec) );

    tree_vec.push_back(1.0, *psi_np1);
    tree_vec.push_back(-1.0, *psi_n);
    FunctionTree<3> *d_psi = add(tree_vec);
    tree_vec.clear();

    double error = sqrt(d_psi->getSquareNorm());
    REQUIRE( error < apply_prec );

    delete V;
    delete Vpsi;
    delete psi_n;
    delete psi_np1;
    delete d_psi;
}

} // namespace
