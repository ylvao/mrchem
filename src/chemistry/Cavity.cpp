#include "Cavity.h"
#include <cmath> 
#include <vector>
#include <array>

namespace mrchem {

Cavity::Cavity(std::vector<std::array<double, 3>> coord, std::vector<double> R, double slope){
  this->pos = coord;
  this->R = R;
  this->d = slope; 
}


double Cavity::evalf(const double *r) const {
  double C = 1.0;
  double s, O;  

  for(int i = 0; i < pos.size(); i++){
    s = std::sqrt(std::pow(pos[i][0] - r[0], 2) + std::pow(pos[i][1] - r[1], 2) + std::pow(pos[i][2] - r[2], 2)) - R[i];
    O = 0.5 * (1 + std::erf(s/d));
    C *= 1 - (1 - O);
  }
  C = 1 - C;
  return C;
}



} //namespace mrchem
