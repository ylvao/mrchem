#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>

#include "Density.h"
#include "Potential.h"
#include "OrbitalSet.h"
#include "DensExp.h"

using namespace std;

/** Compute the density from a Gaussian density expansion
 */
Density::Density(DensExp &densExp) {
    double screening = -log(this->getRelPrec());
    screening = 32.0;
    densExp.calcScreening(screening);
    densExp.setScreen(true);
    this->projectFunction(densExp);
    densExp.setScreen(false);
}

Density::Density(Orbital &orb_i, Orbital &orb_j, int s, bool norm) {
    int occ = orb_i.compareOccupancy(orb_j);
    int spin = orb_j.compareSpin(orb_j);

    if (occ < 0) MSG_ERROR("Mixing occupancies");
    if (spin < 0) MSG_ERROR("Mixing spins");

    if (s == Orbital::Alpha) {
	if (spin == Orbital::Undefined) {
	    occ = 1;
	} else if (spin == Orbital::Beta) {
	    occ = 0;
	}
    } else if (s == Orbital::Beta) {
	if (spin == Orbital::Undefined) {
	    occ = 1;
	} else if (spin == Orbital::Alpha) {
	    occ = 0;
	}
    }
    double norm_i = 1.0;
    double norm_j = 1.0;
    if (norm) {
        norm_i = sqrt(orb_i.getSquareNorm());
	norm_j = sqrt(orb_j.getSquareNorm());
    }
    if (occ > 0) {
        this->mult(1.0*occ/norm_i, orb_i, 1.0/norm_j, orb_j);
    }
}

Density::Density(OrbitalSet &i_orbs, OrbitalSet &j_orbs, int s) {
    int N_i = i_orbs.size();
    int N_j = j_orbs.size();
    if (N_i != N_j) MSG_ERROR("Different number of orbitals");

    vector<FunctionTree<3> *> trees;
    for (int p = 0; p < N_i; p++) {
	Orbital &iOrb_p = i_orbs.getOrbital(p);
	Orbital &jOrb_p = j_orbs.getOrbital(p);
	
	Density *rho_p = new Density(iOrb_p, jOrb_p, s);
	trees.push_back(rho_p);
    }
    this->add(trees);
    while (trees.size() > 0) {
	delete trees.back();
	trees.pop_back();
    }
}

Density::Density(Orbital &phi, Orbital &x, Orbital &y, int s) {
    //Orbital orb_1(phi);
    //Orbital orb_2(phi);

    //orb_1.add(1.0, x, 1.0, phi);
    //orb_2.add(1.0, phi, 1.0, y);

    //orb_1.normalize();
    //orb_2.normalize();

    //Orbital dOrb_1(phi);
    //Orbital dOrb_2(phi);
    //dOrb_1.add(1.0, orb_1, -1.0, phi);
    //dOrb_2.add(-1.0, phi, 1.0, orb_2);

    Density x_phi(x, phi, s, false);
    Density phi_y(phi, y, s, false);
    this->add(1.0, x_phi, 1.0, phi_y);
}

Density::Density(OrbitalSet &orbs, OrbitalSet &x, OrbitalSet &y, int s) {
    vector<double> coefs;
    vector<FunctionTree<3> *> trees;

    int nOrbs = orbs.size();
    if (nOrbs != x.size()) MSG_ERROR("Different number of orbitals");
    if (nOrbs != y.size()) MSG_ERROR("Different number of orbitals");
    for (int p = 0; p < nOrbs; p++) {
        Orbital &phi_p = orbs.getOrbital(p);
        Orbital &x_p = x.getOrbital(p);
        Orbital &y_p = y.getOrbital(p);

	Density *rho_p = new Density(phi_p, x_p, y_p, s);
	trees.push_back(rho_p);
    }

    this->add(trees);
    while (trees.size() > 0) {
	delete trees.back();
	trees.pop_back();
    }
}

/*
        int orbSpin = orb.getSpin();
        int occ = orb.getOccupancy();

	int n = 0;
        if (densSpin == Orbital::Undefined) {
            n = orb.getOccupancy();
        } else if (orbSpin == Orbital::Undefined) {
            n = 1;
        } else if (orbSpin == densSpin) {
            n = 1;
        }

        if (n > 0) {
            orbs.fetch(i);
            Orbital *sqOrb = new Orbital;
            sqOrb->mult(1.0, orb, 1.0, orb);
            orbs.dump(i, false);

            double norm = orb.getSquareNorm();
            if (dump) {
                Density tmpDens(*this);
                this->clear();
                this->add(1.0, tmpDens, 1.0*n/norm, *sqOrb);
                delete sqOrb;
            } else {
                coefs.push_back(1.0*n/norm);
                trees.push_back(sqOrb);
            }
        }
    }
    if (not dump) {
        this->add(coefs, trees);
    }
    for (int n = 0; n < trees.size(); n++) {
        delete trees[n];
        trees[n] = 0;
    }
}
*/

/** Compute spin density from an orbital set
 *
 * Squares each orbital and adds them up weighted by their occupancy
 * \f$ \rho^\sigma(r) = \sum_i occ_i^sigma |\phi_i^\sigma(r)|^2\f$
 * The squared orbitals are automatically normalized to yield the 
 * correct number density, but they are NOT orthogonalized.
 * Undefined spin input computes the full density \f$\alpha + \beta\f$
 * Default spin input is Undefined. If the orbitals are stored on disk 
 * each one is loaded, squared, added and dumped separately. If all 
 * orbitals are in memory they are first all squared and then added 
 * up simultaneously to the density. This makes the process much slower
 * in the case of stored orbitals. 
 */
/*
Density::Density(OrbitalSet &orbs, int densSpin) {
    bool dump = orbs.getDumpToFile();

    vector<double> coefs;
    vector<FunctionTree<3> *> trees;

    for (int i = 0; i < orbs.size(); i++) {
        Orbital &orb = orbs.getOrbital(i);
        int orbSpin = orb.getSpin();
        int occ = orb.getOccupancy();

	int n = 0;
        if (densSpin == Orbital::Undefined) {
            n = orb.getOccupancy();
        } else if (orbSpin == Orbital::Undefined) {
            n = 1;
        } else if (orbSpin == densSpin) {
            n = 1;
        }

        if (n > 0) {
            orbs.fetch(i);
            Orbital *sqOrb = new Orbital;
            sqOrb->mult(1.0, orb, 1.0, orb);
            orbs.dump(i, false);

            double norm = orb.getSquareNorm();
            if (dump) {
                Density tmpDens(*this);
                this->clear();
                this->add(1.0, tmpDens, 1.0*n/norm, *sqOrb);
                delete sqOrb;
            } else {
                coefs.push_back(1.0*n/norm);
                trees.push_back(sqOrb);
            }
        }
    }
    if (not dump) {
        this->add(coefs, trees);
    }
    for (int n = 0; n < trees.size(); n++) {
        delete trees[n];
        trees[n] = 0;
    }
}
*/

/*
Density::Density(OrbitalSet &aOrbs, OrbitalSet &bOrbs, int densSpin) {
    if (aOrbs.size() != bOrbs.size()) {
	MSG_ERROR("Size mismatch");
    }

    bool dump = false;
    if (aOrbs.getDumpToFile() or bOrbs.getDumpToFile()) {
	dump = true;
    }

    vector<double> coefs;
    vector<FunctionTree<3> *> trees;

    for (int i = 0; i < aOrbs.size(); i++) {
        Orbital &aOrb = aOrbs.getOrbital(i);
        Orbital &bOrb = bOrbs.getOrbital(i);

        int orbSpin = aOrb.compareSpin(bOrb);
	if (orbSpin < 0) {
	    MSG_ERROR("Spin mismatch orbital " << i);
	}

        int occ = aOrb.compareOccupancy(bOrb);
	if (occ < 0) {
	    MSG_ERROR("Occupancy mismatch orbital " << i);
	}

	int n = 0;
        if (densSpin == Orbital::Undefined) {
            n = occ;
        } else if (orbSpin == Orbital::Undefined) {
            n = 1;
        } else if (orbSpin == densSpin) {
            n = 1;
        }

        if (n > 0) {
            aOrbs.fetch(i);
            bOrbs.fetch(i);
            Orbital *sqOrb = new Orbital;
            sqOrb->mult(1.0, aOrb, 1.0, bOrb);
            aOrbs.dump(i, false);
            bOrbs.dump(i, false);

            if (dump) {
                Density tmpDens(*this);
                this->clear();
                this->add(1.0, tmpDens, 1.0*n, *sqOrb);
                delete sqOrb;
            } else {
                coefs.push_back(1.0*n);
                trees.push_back(sqOrb);
            }
        }
    }
    if (not dump) {
        this->add(coefs, trees);
    }
    for (int n = 0; n < trees.size(); n++) {
        delete trees[n];
        trees[n] = 0;
    }
}
*/

/** Density copy constructor
 *
 * Makes a complete copy of the density, with tree parameters.
 */
Density::Density(const Density &dens): FunctionTree<3>(dens) {
}

/** Density assignment operator
 *
 * This should only copy the density function itself,
 * not the tree parameters such as precision.
 */
Density &Density::operator=(const Density &dens) {
    if (this == &dens) {
        return *this;
    }
    FunctionTree<3>::operator=(dens);
    return *this;
}

/** Renormalizes the density to yield the given density integral. */
void Density::setNumberDensity(double numDens) {
    double oldCharge = this->integrate();
    *this *= numDens/oldCharge;
}

/** Write the tree structure to disk, for later use.
  * Argument file name will get a ".dens" file extension, and in MPI an
  * additional "-[rank]". */
bool Density::saveTree(const string &file) {
    stringstream fname;
    fname << file;
    if (this->isScattered()) {
        fname << "-" << this->getRankId();
    }
    fname << ".dens";
    ofstream ofs(fname.str().c_str(), ios_base::binary);
    if (ofs == 0) {
        MSG_FATAL("Could not open file for writing: " << file);
    }
    boost::archive::binary_oarchive oa(ofs);
    this->purgeGenNodes();
    oa << *this;
    return true;
}

/** Read a previously stored tree structure from disk.
  * Argument file name will get a ".dens" file extension, and in MPI an
  * additional "-[rank]". */
bool Density::loadTree(const string &file) {
    stringstream fname;
     
    fname << file;
    if (node_group.size() > 1 and this->isBuildDistributed()) {
        fname << "-" << this->getRankId();
    }
    fname << ".dens";
    ifstream ifs(fname.str().c_str(), ios_base::binary);
    if (not ifs) {
        return false;
    }
    boost::archive::binary_iarchive ia(ifs);
    ia >> *this;
    return true;
}
