#include "logsumexp.h"

double logsumexp(const Eigen::VectorXd &x){
    double max = x.maxCoeff();
    return max + std::log((x.array() - max).exp().sum());
}