#include <Eigen/Dense>
#include "PolyInterpolator.h"

double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x){
    double xm2 = x_in(0);
    double xm1 = x_in(1);
    double x00 = x_in(2);
    double xp1 = x_in(3);
    double xp2 = x_in(4);

    double ym2 = y_in(0);
    double ym1 = y_in(1);
    double y00 = y_in(2);
    double yp1 = y_in(3);
    double yp2 = y_in(4);
    return ym2 + (x - xm2)*((ym1 - ym2)/(xm1 - xm2) +
        (x - xm1)*(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2) +
        (x - x00)*((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/(-xm2 + xp1) +
        ((x - xp1)*(-((-(((y00 - ym1)/(x00 - xm1) + (-ym1 + ym2)/(xm1 - xm2))/(x00 - xm2)) +
        ((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/(-xm1 + xp1))/
        (-xm2 + xp1)) + (-(((-y00 + ym1)/(x00 - xm1) + (y00 - yp1)/(x00 - xp1))/
        (-xm1 + xp1)) + ((-y00 + yp1)/(x00 - xp1) + (yp1 - yp2)/(xp1 - xp2))/(-x00 + xp2))/
        (-xm1 + xp2)))/(-xm2 + xp2))));
}