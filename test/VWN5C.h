/** 
 *
 * \date Jul 13, 2010
 * \author Shubham Nayyar \n
 *		   CTCC, University of Tromsø
 *
 *
 */

#ifndef VWN5C_H
#define VWN5C_H

#include "math_consts.h"
#include "LDAFunc.h"
#include "ProjectedNode.h"
#include "QuadratureCache.h"

namespace vwn {
	static double para[] = {-0.10498, 0.0621814, 3.72744, 12.9352};
//	static double ferro[] = {-0.325, 0.0310907, 7.06042, 18.0578};
//	static double inter[] = {-0.0047584, -pow(3*M_PI*M_PI,-1.0), 1.13107, 13.0045};

//	ulfek: second elements are multiplied by 2 wrt molpro manual
//	ulfek: in Dalton and Dirac inter[1] has too few decimals, should  be -1/(3pi^2)

//	double  vwn_a(double p[]) {
//		return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
//	}
//
//	double vwn_b(double p[]) {
//		return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
//	}
//
//	double vwn_c(double p[]) {
//		return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
//				((p[0]*p[0] + p[0]*p[2] + p[3])*
//				sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
//	}


//	above to be calculated only when p(alpha)!=p(beta)
	const double vwn_a = para[0]*para[2]/(para[0]*para[0] +
						 para[0]*para[2] + para[3]) - 1;
	const double vwn_b = 2*(para[0]*para[2]/(para[0]*para[0] +
						 para[0]*para[2] + para[3]) - 1) + 2;
	const double vwn_c = 2*para[2]*(1/sqrt(4*para[3] - para[2]*para[2]) -
						 para[0]/((para[0]*para[0] + para[0]*para[2] + para[3])*
						 sqrt(4*para[3] - para[2]*para[2])/(para[2] + 2*para[0])));

	static double vwn_x(const double &s, double p[]) {
		return s*s + p[2]*s + p[3];
	}

	static double vwn_y(const double &s, const double p[]) {
		return s - p[0];
	}

	static double vwn_z(const double &s, double p[]) {
		return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
	}

	static double vwn_f(const double &s, double p[]) {
		return 0.5*p[1]*(2*log(s) + vwn_a*log(vwn_x(s, p)) -
			vwn_b*log(vwn_y(s, p)) +
			vwn_c*atan(vwn_z(s, p)));
	}

//	zeta ,g and others are not calculated as p(alpha)-p(beta) = 0

	static double VWN5C (double  d) {
		if (d < MachineZero) {
			return 0.0;
		}
//		using namespace vwn;
//		double zeta = 0;
		double s = 0.25*pow((3*pow(4.0,5)*(1/(M_PI *fabs(d)))),(1.0/6));
//		Constant is (2^1/3-1)^-1/2
//		double g = 1.92366105093154*(pow(1 + zeta,4.0/3) + pow(1 - zeta,4.0/3) - 2);
//		double dd = g*((vwn_f(s, ferro) - vwn_f(s, para))*zeta +
//		vwn_f(s, inter)*(1 - zeta)*(9.0/4.0*(pow(2,1.0/3.0)-1)));
//		Dalton (for some reason??) uses a value:    0.584822305543806
//		The real value should be 9/4 (2^(1/3) -1) = 0.5848223622634646
		return d*vwn_f(s, para);
	}
}

class VWN5C : public LDAFunc<3> {
public:
	VWN5C() { }
	~VWN5C() { }

protected:
	void transformValues(Eigen:: VectorXd &A) {
		for (int i = 0; i < A.size(); i++) {
			A[i] = vwn::VWN5C(A[i]);
		}
	}
};

#endif
