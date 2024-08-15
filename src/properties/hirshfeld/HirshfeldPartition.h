#pragma once

#include <Eigen/Dense>
#include <mrchem.h>
#include "chemistry/Molecule.h"
#include <string>

class HirshfeldPartition{

    public:
        HirshfeldPartition(const mrchem::Molecule &mol, std::string data_dir);

        mrcpp::ComplexFunction getHirshfeldPartitionFunction(int index, double prec) const;

    protected:

};