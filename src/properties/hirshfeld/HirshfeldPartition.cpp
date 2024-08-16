#include "HirshfeldPartition.h"
#include "properties/hirshfeld/logsumexp.h"

HirshfeldPartition::HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir) {

    this->nucs = std::make_shared<mrchem::Nuclei>(mol.getNuclei());
    this->nNucs = this->nucs->size();

}



mrcpp::ComplexFunction HirshfeldPartition::getHirshfeldPartitionFunction(int index, double prec) const {
    
}

double HirshfeldPartition::lseLogDens(mrcpp::Coord<3> &r) const {
    Eigen::VectorXd lseLogDens_r(this->nNucs);
    double rr;
    mrcpp::Coord<3> nucPos;
    for (int i = 0; i < this->nNucs; i++) {
        nucPos = this->nucs->at(i).getCoord();
        rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
            + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
            + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
        lseLogDens_r(i) = this->logDensities[i].evalf(rr);
    }
    return logsumexp(lseLogDens_r);
}

double HirshfeldPartition::evalf(mrcpp::Coord<3> &r, int iAt) const {
    mrcpp::Coord<3> nucPos = this->nucs->at(iAt).getCoord();
    double rr = std::sqrt((r[0] - nucPos[0]) * (r[0] - nucPos[0])
        + (r[1] - nucPos[1]) * (r[1] - nucPos[1])
        + (r[2] - nucPos[2]) * (r[2] - nucPos[2]));
    return std::exp(this->logDensities[iAt].evalf(rr) - this->lseLogDens(r));
}