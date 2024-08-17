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
    HirshfeldRadInterpolater(const std::string element, std::string data_dir, bool writeToFile = false);

    /**
     * @brief Evaluate the interpolated function at a given point
     * @param r The point at which to evaluate
     * @return The interpolated value at the given point
    */
    double evalf(const double &r) const;

    /**
     * @brief Integrates the atomic charge density to get the charge of the tabulated 
     * atomic density. The numeric is performed from 0 to 20 bohr. Only useful for debugging.
     */
    double getNorm() const;

protected:
    /**
     * @brief The interpolator for the atomic density
     */
    std::shared_ptr<PolyInterpolator> lnRho = nullptr;
    /**
     * @brief Write the interpolated density to a file for debugging.
     * @param path The path to the file to write the interpolated density to.
     */
    void writeInterpolatedDensity(const std::string path);
};
