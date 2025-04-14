#pragma once
#include <MRCPP/MWOperators>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <xc.h>  // The addition of LibXC header
#include "MRDFT.h"

#include "LibXC.h"
#include "LibXC.cpp"

namespace mrdft {

// Structure to hold LibXC functional data
struct LibXCData {
    int func_id;             // Functional ID in LibXC
    double weight;           // Weight for this functional
    xc_func_type func;       // LibXC function type
    bool initialized;        // If the functional is initialized
};

class Factory final {
public:
    Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA);
    ~Factory();
    
    void setSpin(bool s) { spin = s; }
    void setOrder(int k) { order = k; }
    void setUseGamma(bool g) { gamma = g; }
    void setLogGradient(bool lg) { log_grad = lg; }
    void setDensityCutoff(double c) { cutoff = c; }
    void setDerivative(const std::string &n) { diff_s = n; }
    
    // New methods for LibXC
    void setFunctional(const std::string &name, double weight = 1.0);
    bool isGGA() const;
    bool isHybrid() const;
    double getHybridCoeff() const;
    
    std::unique_ptr<MRDFT> build();
    
private:
    int order{1};
    bool spin{false};
    bool gamma{false};
    bool log_grad{false};
    double cutoff{-1.0};
    std::string diff_s{"abgv_00"};
    const mrcpp::MultiResolutionAnalysis<3> mra;
    
    // LibXC data
    std::vector<LibXCData> functionals;
    
    // Helper methods
    int mapFunctionalName(const std::string &name) const;
    void cleanupFunctionals();
    
    std::unique_ptr<mrcpp::DerivativeOperator<3>> diff_p;
};

} // namespace mrdft