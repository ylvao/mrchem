#include "HirshfeldInterpolator.h"
#include <filesystem>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

// Function to read atomic density data from a file
void readAtomicDensity(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &rhoGrid) {
    std::vector<double> r, rho;
    std::ifstream file(path);
    std::string line;
    double r_, rho_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> r_ >> rho_;
        if (rho_ > 0) {
            r.push_back(r_);
            rho.push_back(rho_);
        }
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    rhoGrid = Eigen::Map<Eigen::VectorXd>(rho.data(), rho.size());
}

// Constructor
HirshfeldRadInterpolater::HirshfeldRadInterpolater(const std::string element, std::string data_dir) {
    Eigen::VectorXd rGrid;
    Eigen::VectorXd rhoGrid;

    std::string filename = data_dir + '/' + element + ".density";

    readAtomicDensity(filename, rGrid, rhoGrid);

    rhoGrid = rhoGrid.array().log() / (4 * M_PI); // todo multiply by 4pi*r^2 or some other magic number

    lnRho = std::make_shared<PolyInterpolator>(rGrid, rhoGrid);
}

// Function to evaluate the interpolated function
double HirshfeldRadInterpolater::evalf(const double &r) const {
    double y;
    return lnRho->evalf(r);
}
