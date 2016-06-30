#ifndef DIPOLEMOMENT_H
#define DIPOLEMOMENT_H

#include <Eigen/Core>

#include "TelePrompter.h"

class Nuclei;
class OrbitalVector;
class DipoleOperator;

class DipoleMoment {
public:
    DipoleMoment();
    virtual ~DipoleMoment();

    Eigen::VectorXd get() const { return this->dipole_nuc + this->dipole_el; }
    Eigen::VectorXd getNuclear() const { return this->dipole_nuc; }
    Eigen::VectorXd getElectronic() const { return this->dipole_el; }

    void compute(int d, DipoleOperator &mu, const Nuclei &nuc);
    void compute(int d, DipoleOperator &mu, OrbitalVector &phi);

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
