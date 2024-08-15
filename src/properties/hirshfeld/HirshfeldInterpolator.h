#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <filesystem>
#include "PolyInterpolator.h"

void readAtomicDensity(const std::string path, Eigen::VectorXd &rGrid, Eigen::VectorXd &rhoGrid){
    std::vector<double> r, rho;
    std::ifstream file(path);
    std::string line;
    double r_, rho_;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        iss >> r_ >> rho_;
        r.push_back(r_);
        rho.push_back(rho_);
    }
    file.close();
    rGrid = Eigen::Map<Eigen::VectorXd>(r.data(), r.size());
    rhoGrid = Eigen::Map<Eigen::VectorXd>(rho.data(), rho.size());
}

class HirshfeldRadInterpolater {

    public:
    /**
     * @brief Construct a new Rad Interpolater object
     * @param element The element for which the ZORA potential is to be interpolated
     * @param data_dir The directory containing the ZORA potential data
    */
    RadInterpolater(const std::string element, std::string data_dir){
        Eigen::VectorXd rGrid;
        Eigen::VectorXd rhoGrid;


        this->mode = mode;
        std::string filename = data_dir + '/' + element + ".dat";

        readAtomicDensity(filename, rGrid, rhoGrid);

        rhoGrid = Eigen::log(rhoGrid); // todo multiply by 4pi*r^2 or some other magic number

        lnRho = std::make_shared<PolyInterpolator>(rGrid, rhoGrid);

    }

    double evalf(const double &r) const {
        double y, yp;
        return lnRho->eval(r, y);
    }

    protected:
    std::shared_ptr<PolyInterpolator> lnRho = nullptr;

};
