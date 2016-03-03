/** 
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#ifndef DIRACEXCHANGED1_H
#define DIRACEXCHANGED1_H

#include "math_consts.h"
#include "LDAFunc.h"
#include "ProjectedNode.h"
#include "QuadratureCache.h"

class DiracExchangeD1 : public LDAFunc<3> {
public:
	DiracExchangeD1() { }
	~DiracExchangeD1() { }
protected:
	void transformValues(Eigen::VectorXd &val) {
		for (int i = 0; i < val.size(); i++) {
			val[i] = cbrt(val[i]);
		}
		val *= 4.0/3.0 * C_x;
	}
};

#endif // DIRACEXCHANGED1_H
