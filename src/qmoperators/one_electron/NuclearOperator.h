#pragma once

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class NuclearPotential final : public QMPotential {
public:
    NuclearPotential(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false);
    ~NuclearPotential() override { free(NUMBER::Total); }

    void setup(double prec) override { setApplyPrec(prec); }
    void clear() override { clearApplyPrec(); }

    Nuclei &getNuclei() { return this->func.getNuclei(); }
    const Nuclei &getNuclei() const { return this->func.getNuclei(); }
    double evalf(const mrcpp::Coord<3> &r) { return this->func.evalf(r); }

private:
    NuclearFunction func;

    void allreducePotential(double prec, QMFunction &V_loc);
};

class NuclearOperator final : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false) {
        r_m1 = std::make_shared<NuclearPotential>(nucs, proj_prec, smooth_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
    }

    Nuclei &getNuclei() { return this->r_m1->getNuclei(); }
    const Nuclei &getNuclei() const { return this->r_m1->getNuclei(); }

    double trace(const Nuclei &nucs);

    using RankZeroTensorOperator::trace;

private:
    std::shared_ptr<NuclearPotential> r_m1{nullptr};
};

} // namespace mrchem
