#include <vector>
#include <array>
#include "MRCPP/MWFunctions"

namespace mrchem {


class Cavity : public mrcpp::RepresentableFunction<3> {
public: 
  Cavity(std::vector<std::vector<double>> coord, std::vector<double> R, double slope);
  double evalf(const double *r) const;
  double evalf(const std::array<double, 3> &r) const {return 0;}
protected:
  std::vector<std::vector<double>> pos;
  std::vector<double> R;
  double d;
};

}

