#include "DipoleMoment.h"
#include "OrbitalVector.h"
#include "DipoleOperator.h"

DipoleMoment::DipoleMoment() {
    this->dipole_nuc.setZero();
    this->dipole_el.setZero();
}

DipoleMoment::~DipoleMoment() {
}

void DipoleMoment::compute(int d, DipoleOperator &mu, const Nuclei &nucs) {
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc_k = nucs[k];
        double Z = nuc_k.getCharge();
        this->dipole_nuc(d) += Z*mu(nuc_k);
    }
}

void DipoleMoment::compute(int d, DipoleOperator &mu, OrbitalVector &orbs) {
    for (int n = 0; n < orbs.size(); n++) {
        Orbital &orb_n = orbs.getOrbital(n);
        double occ = (double) orb_n.getOccupancy();
        this->dipole_el(d) -= occ*mu(orb_n, orb_n);
    }
}
