#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "XCPotential.h"
#include "XCPotentialD2.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/density_utils.h"

/** @brief Return FunctionTree for the input density from the XCFunctional
 *
 * @param[in] type Which density to return (alpha, beta or total)
 */
FunctionTree<3> &XCPotential::getDensity(int spin) {
    if (spin == DENSITY::Total) return this->functional->getDensity(mrdft::DensityType::Total);
    if (spin == DENSITY::Alpha) return this->functional->getDensity(mrdft::DensityType::Alpha);
    if (spin == DENSITY::Beta)  return this->functional->getDensity(mrdft::DensityType::Beta);
    MSG_FATAL("Invalid density type");
}
