#pragma once
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

double polynomialInterpolate5_deriv(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

class PolyInterpolator {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    int n;
    double xmin;
    double xmax;
    double ypxmax;
    public:
    PolyInterpolator(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in){
        x = x_in;
        y = y_in;
        n = x_in.size();
        xmin = x_in(0);
        xmax = x_in(x_in.size() - 1);

        Eigen::VectorXd x_in_poly(5), y_in_poly(5);
        for (int i = n - 5; i < n; i++) {
            x_in_poly(i) = x(i);
            y_in_poly(i) = y(i);
        }
        ypxmax = polynomialInterpolate5_deriv(x_in_poly, y_in_poly, xmax);
    }
    
    double evalf(const double &x) const {
        double y;
        if (x > xmax) { // linear extrapolation
                y = this->y(n - 1) + ypxmax*(x - xmax);
            return y;
        }

        // Find the interval in which x lies using binary search
        int i = 0;
        int j = this->x.size() - 1;
        while (j - i > 1) {
            int k = (i + j)/2;
            if (x < this->x(k)) {
                j = k;
            } else {
                i = k;
            }
        }
        Eigen::VectorXd x_in(5), y_in(5);
        if (i == 0) i = 2;
        if (i == 1) i = 2;
        if (i == this->x.size() - 1) i = this->x.size() - 3;
        if (i == this->x.size() - 2) i = this->x.size() - 3;
        for (int k = 0; k < 5; k++) {
            x_in(k) = this->x(i - 2 + k);
            y_in(k) = this->y(i - 2 + k);
        }
        y = polynomialInterpolate5(x_in, y_in, x);
        return y;
    }
};