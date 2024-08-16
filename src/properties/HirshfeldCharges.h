#pragma once

#include <nlohmann/json.hpp>

#include "mrchem.h"

#include "utils/math_utils.h"
#include "utils/print_utils.h"

namespace mrchem {

class HirshfeldCharges final {
public:

    DoubleVector getVector() const { return this->hirshfeld_charges; }

    void setVector(const DoubleVector &v) { this->hirshfeld_charges = v; }

    void print(const std::string &id) const {
        mrcpp::print::header(0, "Hirshfeld Charges (" + id + ")");
        print_utils::vector(0, "Hirshfeld Charges", getVector(), 1e-12);
        mrcpp::print::separator(0, '-');
        print_utils::scalar(0, "Sum of Hirshfeld charges", getVector().sum(), "(au)", 1e-12);
    }

    nlohmann::json json() const {
        return {{"total", print_utils::eigen_to_vector(getVector(), 1.0e-12)}};
    }

protected:
    DoubleVector hirshfeld_charges;
};

} // namespace mrchem