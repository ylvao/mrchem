#pragma once

#include <Eigen/Dense>
#include <mrchem.h>
#include "chemistry/Molecule.h"
#include <string>
#include "chemistry/Nucleus.h"
#include "properties/hirshfeld/HirshfeldInterpolator.h"

class HirshfeldPartition{

    public:
        HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir);

        mrcpp::ComplexFunction getHirshfeldPartitionFunction(int index, double prec) const;

    protected:

    double evalf(const mrcpp::Coord<3> &r, int iAt) const;
    double lseLogDens(const mrcpp::Coord<3> &r) const;

    std::shared_ptr<mrchem::Nuclei> nucs{nullptr};

    int nNucs;

    std::vector<HirshfeldRadInterpolater> logDensities;


};