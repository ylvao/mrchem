/** 
 *
 * \date Jul 13, 2010
 * \author Stig Rune Jensen \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#ifndef DIRACEXCHANGE_H
#define DIRACEXCHANGE_H

#include "math_consts.h"
#include "LDAFunc.h"
#include "ProjectedNode.h"
#include "QuadratureCache.h"

class DiracExchange : public LDAFunc<3> {
public:
	DiracExchange() { }
	~DiracExchange() { }
protected:
	void transformValues(Eigen::VectorXd &val) {
		for (int i = 0; i < val.size(); i++) {
            if (val[i] > 1.0e-12) {
                val[i] = pow(val[i], 4.0/3.0);
            }
            else {
                val[i] = 0.0;
            }
		}
        val *= C_x;
	}
};

#endif // DIRACEXCHANGE_H
