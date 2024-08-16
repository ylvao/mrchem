#pragma once

#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>
#include <memory>
#include "PolyInterpolator.h"

/**
 * @brief Reads atomic density data from a file
 * @param path The path to the file containing the data
 * @param rGrid The grid for the radial distances
 * @param rhoGrid The grid for the atomic density values
 */
void readAtomicDensity(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &rhoGrid);

class HirshfeldRadInterpolater {

public:
    /**
     * @brief Construct a new Rad Interpolater object
     * @param element The element for which the ZORA potential is to be interpolated
     * @param data_dir The directory containing the ZORA potential data
    */
    HirshfeldRadInterpolater(const std::string element, std::string data_dir);

    /**
     * @brief Evaluate the interpolated function at a given point
     * @param r The point at which to evaluate
     * @return The interpolated value at the given point
    */
    double evalf(const double &r) const;

protected:
    std::shared_ptr<PolyInterpolator> lnRho = nullptr;
    void writeInterpolatedDensity(const std::string path);
};
