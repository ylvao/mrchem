#pragma once
#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <filesystem>
#include <iostream>

/**
 * @brief Interpolate a 5th order polynomial through 5 points and evaluate at x.
 * @param x_in x values of the 5 points
 * @param y_in y values of the 5 points
 * @param x x value at which to evaluate the polynomial
 */
double polynomialInterpolate5(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

/**
 * @brief Evaluate the derivative of a 5th order polynomial through 5 points at x.
 * @param x_in x values of the 5 points
 * @param y_in y values of the 5 points
 * @param x x value at which to evaluate the derivative of the polynomial
 */
double polynomialInterpolate5_deriv(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in, double x);

/**
 * @brief Class to interpolate a functions using a 5th order polynomials.
 */
class PolyInterpolator {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    int n;
    double xmin;
    double xmax;
    double ypxmax;
    public:
    /**
     * @brief Constructor
     * @param x_in x values of the points to interpolate
     * @param y_in y values of the points to interpolate
     */
    PolyInterpolator(Eigen::VectorXd &x_in, Eigen::VectorXd &y_in){
        x = x_in;
        y = y_in;
        n = x_in.size();
        xmin = x_in(0);
        xmax = x_in(x_in.size() - 1);

        Eigen::VectorXd x_in_poly(5), y_in_poly(5);
        int j = 0;
        for (int i = n - 5; i < n; i++) {
            x_in_poly(j) = x(i);
            y_in_poly(j) = y(i);
            j++;
        }
        ypxmax = polynomialInterpolate5_deriv(x_in_poly, y_in_poly, xmax);
    }
    
    /**
     * @brief Evaluate the interpolated function at x. No extrapolation for x < xmin, linear extrapolation for x > xmax.
     * @param x x value at which to evaluate the function
     */
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