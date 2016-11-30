#ifndef DIPOLEMOMENT_H
#define DIPOLEMOMENT_H

#include <Eigen/Core>

#include "TelePrompter.h"
#include "OrbitalVector.h"
#include "DipoleOperator.h"

class DipoleMoment {
public:
    DipoleMoment() {
        this->dipole_nuc.setZero();
        this->dipole_el.setZero();
    }
    virtual ~DipoleMoment() { }

    Eigen::Vector3d get() const { return this->dipole_nuc + this->dipole_el; }
    Eigen::Vector3d getNuclear() const { return this->dipole_nuc; }
    Eigen::Vector3d getElectronic() const { return this->dipole_el; }

    void compute(int d, DipoleOperator &mu, const Nuclei &nucs) {
        for (int k = 0; k < nucs.size(); k++) {
            const Nucleus &nuc_k = nucs[k];
            double Z = nuc_k.getCharge();
            this->dipole_nuc(d) += Z*mu(nuc_k);
        }
    }
    void compute(int d, DipoleOperator &mu, OrbitalVector &orbs) {
        for (int n = 0; n < orbs.size(); n++) {
            Orbital &orb_n = orbs.getOrbital(n);
            double occ = (double) orb_n.getOccupancy();
            this->dipole_el(d) -= occ*mu(orb_n, orb_n);
        }
    }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dipole) {
        int oldPrec = TelePrompter::setPrecision(10);
        double au = dipole.get().norm();
        double debye = au/0.393430307;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                         Dipole moment                      "<<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<<" Length of vector:   (au)    " << std::setw(30) << au        <<std::endl;
        o<<"                     (Debye) " << std::setw(30) << debye     <<std::endl;
        o<<"------------------------------------------------------------"<<std::endl;
        o<< std::setw(19) << dipole.get()(0);
        o<< std::setw(20) << dipole.get()(1);
        o<< std::setw(20) << dipole.get()(2) << std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    Eigen::Vector3d dipole_nuc;
    Eigen::Vector3d dipole_el;
};

#endif // DIPOLEMOMENT_H
