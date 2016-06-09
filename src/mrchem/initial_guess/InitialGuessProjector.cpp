#include "InitialGuessProjector.h"
#include "GaussExp.h"
#include "Intgrl.h"
#include "OrbitalExp.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "TelePrompter.h"
#include "Timer.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

//void OrbitalVector::readOrbitals(const OrbitalVector &orbs) {
//    boost::timer rolex;
//    rolex.restart();
//    int oldPrec = TelePrompter::setPrecision(15);
//    printout(0, "\n\n================ Setting up starting guess ");
//    printout(0, "=================\n\n");

//    for (int i = 0; i < this->size(); i++) {
//    Orbital *thisOrb = this->getOrbitalPtr(i);
//        if (thisOrb == 0) MSG_ERROR("Orbital not initialized");
//    const Orbital &thatOrb = orbs.getOrbital(i);
//    *thisOrb = thatOrb;
//        printout(0, "Orbital " << setw(3) << i);
//        println(0, " squareNorm: " << setw(36) << thisOrb->getSquareNorm());
//    }
//    TelePrompter::setPrecision(5);
//    printout(0, "\n================ Elapsed time: ");
//    println(0, rolex.elapsed() << " =================\n");
//    TelePrompter::setPrecision(oldPrec);
//}

void InitialGuessProjector::operator()(OrbitalVector &orbs,
                                       const string &bf,
                                       const string &mo) {
    Timer rolex;
    rolex.restart();
    TelePrompter::printHeader("Setting up occupied orbitals");
    println(0, "    n  Spin  Occ                           SquareNorm");
    TelePrompter::printSeparator('-');

    OrbitalExp *moExp = readOrbitalExpansion(bf, mo);
    for (int i = 0; i < orbs.size(); i++) {
        Orbital &mwOrb = orbs.getOrbital(i);
        GaussExp<3> &gtOrb = moExp->getOrbital(i);
        mwOrb.clear();
        mwOrb.real = G();
        Q(*mwOrb.real, gtOrb);
        printout(0, setw(5) << i);
        printout(0, setw(5) << mwOrb.printSpin());
        printout(0, setw(5) << mwOrb.getOccupancy());
        printout(0, setw(44) << mwOrb.getSquareNorm() << endl);
    }
    delete moExp;

    TelePrompter::printFooter(rolex, 2);
}

void InitialGuessProjector::operator()(OrbitalVector &orbs,
                                       const string &bf,
                                       const string &mo_a,
                                       const string &mo_b) {
    NOT_IMPLEMENTED_ABORT;
//    boost::timer rolex;
//    rolex.restart();
//    int oldPrec = TelePrompter::setPrecision(15);
//    printout(0, "\n\n================ Setting up starting guess ");
//    printout(0, "=================\n\n");

//    int n_a = 0;
//    int n_b = 0;
//    OrbitalExp *moExp_a = readOrbitalExpansion(bf, mo_a);
//    OrbitalExp *moExp_b = readOrbitalExpansion(bf, mo_b);
//    for (int i = 0; i < this->size(); i++) {
//        Orbital *mwOrb = this->getOrbitalPtr(i);
//        if (mwOrb == 0) MSG_ERROR("Orbital not initialized");
//        GaussExp<3> *gtOrb = 0;
//        if (mwOrb->getSpin() == Orbital::Alpha) {
//            gtOrb = &moExp_a->getOrbital(n_a);
//            n_a++;
//        }
//        if (mwOrb->getSpin() == Orbital::Beta) {
//            gtOrb = &moExp_b->getOrbital(n_b);
//            n_b++;
//        }
//        if (mwOrb->getSpin() == Orbital::Paired) {
//            NOT_IMPLEMENTED_ABORT;
//        }
//        //mwOrb.setRelPrec(1.0e-3);
//        mwOrb->projectFunction(*gtOrb);
//        printout(0, "Orbital " << setw(3) << i);
//        println(0, " squareNorm: " << setw(36) << mwOrb->getSquareNorm());
//        //mwOrb.normalize();
//    }
//    delete moExp_a;
//    delete moExp_b;

//    TelePrompter::setPrecision(5);
//    printout(0, "\n================ Elapsed time: ");
//    println(0, rolex.elapsed() << " =================\n");
//    TelePrompter::setPrecision(oldPrec);
}

//void OrbitalVector::readVirtuals(const string &bf, const string &mo, int n_occ) {
//    boost::timer rolex;
//    rolex.restart();
//    int oldPrec = TelePrompter::setPrecision(15);
//    printout(0, "\n\n=============== Setting up virtual orbitals ");
//    printout(0, "================\n\n");

//    OrbitalExp *moExp = readOrbitalExpansion(bf, mo);
//    for (int a = n_occ; a < moExp->size(); a++) {
//    GaussExp<3> &gtOrb = moExp->getOrbital(a);
//        Orbital *orb_a = new Orbital(2, Orbital::Paired);
//    orb_a->projectFunction(gtOrb);
//        printout(0, "Orbital " << setw(3) << a);
//        println(0, " squareNorm: " << setw(36) << orb_a->getSquareNorm());
//        this->orbitals.push_back(orb_a);
//    }
//    delete moExp;
//    TelePrompter::setPrecision(5);
//    printout(0, "\n================ Elapsed time: ");
//    println(0, rolex.elapsed() << " =================\n");
//    TelePrompter::setPrecision(oldPrec);
//}


OrbitalExp* InitialGuessProjector::readOrbitalExpansion(const string &bf, const string &mo) {
    MatrixXd MO = MathUtils::readMatrixFile(mo);
    MatrixXd MO_T = MO.transpose();

    Intgrl intgrl(bf);
    OrbitalExp *moExp = new OrbitalExp(intgrl);
    moExp->rotate(MO_T);
    return moExp;
}

