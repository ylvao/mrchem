#include <vector>
#include <array>
#include "MRCPP/MWFunctions"

using namespace mrcpp;

namespace mrchem {


class Cavity : public mrcpp::RepresentableFunction<3> {
public: 
  Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope);
  double evalf(const double *r) const override; 
  double evalf(const std::array<double, 3> &r) const {return evalf(r.data());};
protected:
  std::vector<std::array<double, 3>> pos;
  std::vector<double> R;
  double d;
};

}

