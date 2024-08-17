#pragma once
#include <Eigen/Dense>

/**
 * @brief Compute the log of the sum of exponentials of the elements of a vector.
 * Numerically stable for large and small x.
 * @param x Vector of values where the logsumexp is to be computed.
 */
double logsumexp(const Eigen::VectorXd &x);