#ifndef SPINSPINCOUPLING_H
#define SPINSPINCOUPLING_H

#include <Eigen/Core>

#include "Nucleus.h"

class SpinSpinCoupling {
public:
    SpinSpinCoupling(const Nucleus &n_k, const Nucleus &n_l)
        : nuc_K(n_k), nuc_L(n_l) {
        this->diamagnetic.setZero();
        this->paramagnetic.setZero();
    }

    virtual ~SpinSpinCoupling() { }

    const Nucleus &getNucleusK() const { return this->nuc_K; }
    const Nucleus &getNucleusL() const { return this->nuc_L; }

    Eigen::Matrix3d get() const { return this->diamagnetic + this->paramagnetic; }
    Eigen::Matrix3d& getDiamagnetic() { return this->diamagnetic; }
    Eigen::Matrix3d& getParamagnetic() { return this->paramagnetic; }

    friend std::ostream& operator<<(std::ostream &o, const SpinSpinCoupling &sscc) {
        double isoDShz = sscc.diamagnetic.trace()/3.0;
        double isoPShz = sscc.paramagnetic.trace()/3.0;
        double isoTShz = isoDShz + isoPShz;

        int oldPrec = TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"================= Spin-Spin Coupling tensor ================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(5);
        o<<std::setw(3)  << sscc.getNucleusK().getElement().getSymbol();
        o<<std::setw(26) << sscc.getNucleusK().getCoord()[0];
        o<<std::setw(15) << sscc.getNucleusK().getCoord()[1];
        o<<std::setw(15) << sscc.getNucleusK().getCoord()[2];
        o<<std::endl;
        o<<std::setw(3)  << sscc.getNucleusL().getElement().getSymbol();
        o<<std::setw(26) << sscc.getNucleusL().getCoord()[0];
        o<<std::setw(15) << sscc.getNucleusL().getCoord()[1];
        o<<std::setw(15) << sscc.getNucleusL().getCoord()[2];
        o<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(10);
        o<<"                                                            "<<std::endl;
        o<<"-------------------- Isotropic averages --------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<<" Total            (Hz)        " << std::setw(30) << isoTShz  <<std::endl;
        o<<" Diamagnetic      (Hz)        " << std::setw(30) << isoDShz  <<std::endl;
        o<<" Paramagnetic     (Hz)        " << std::setw(30) << isoPShz  <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"-------------------------- Total ---------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< sscc.diamagnetic + sscc.paramagnetic                            <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Diamagnetic ------------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< sscc.diamagnetic                                               <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"----------------------- Paramagnetic -----------------------"<<std::endl;
        o<<"                                                            "<<std::endl;
        o<< sscc.paramagnetic                                              <<std::endl;
        o<<"                                                            "<<std::endl;
        o<<"============================================================"<<std::endl;
        o<<"                                                            "<<std::endl;
        TelePrompter::setPrecision(oldPrec);
        return o;
    }
protected:
    const Nucleus nuc_K;
    const Nucleus nuc_L;
    Eigen::Matrix3d diamagnetic;
    Eigen::Matrix3d paramagnetic;
};

#endif // SPINSPINCOUPLING_H
