#ifndef NMRSHIELDING_H
#define NMRSHIELDING_H

#include <Eigen/Core>

#include "Nucleus.h"

class NMRShielding {
public:
    NMRShielding(const Nucleus &n) : nuc(n) {
        this->diamagnetic.setZero();
        this->paramagnetic.setZero();
    }
    virtual ~NMRShielding() { }

    const Nucleus& getNucleus() const { return this->nuc; }

    Eigen::Matrix3d get() const { return this->diamagnetic + this->paramagnetic; }
    Eigen::Matrix3d& getDiamagnetic() { return this->diamagnetic; }
    Eigen::Matrix3d& getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const NMRShielding &nmr) {
        double isoDSppm = nmr.diamagnetic.trace()/3.0;
        double isoPSppm = nmr.paramagnetic.trace()/3.0;
        double isoTSppm = isoDSppm + isoPSppm;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"=================== NMR shielding tensor ==================="<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(5);
        o<<std::setw(3)  << nmr.getNucleus().getElement().getSymbol();
        o<<std::setw(26) << nmr.getNucleus().getCoord()[0];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[1];
        o<<std::setw(15) << nmr.getNucleus().getCoord()[2];
        o<<std::endl;
        TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (ppm)       " << std::setw(30) << isoTSppm <<std::endl;
        o<<" Diamagnetic      (ppm)       " << std::setw(30) << isoDSppm <<std::endl;
        o<<" Paramagnetic     (ppm)       " << std::setw(30) << isoPSppm <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Total ---------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< nmr.diamagnetic + nmr.paramagnetic                            <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Diamagnetic ------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< nmr.diamagnetic                                               <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Paramagnetic -----------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< nmr.paramagnetic                                              <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc;
    Eigen::Matrix3d diamagnetic;
    Eigen::Matrix3d paramagnetic;
};

#endif // NMRSHIELDING_H
