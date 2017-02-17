#ifndef H_BM_DIA_H
#define H_BM_DIA_H

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "CubicPotential.h"

class H_BM_dia : public RankTwoTensorOperator<3,3> {
public:
    H_BM_dia(const double *o = 0, const double *k = 0)
            : o_x(0, o), o_y(1, o), o_z(2, o),
              k_x(0, k), k_y(1, k), k_z(2, k),
              r_m3(1.0, k) {
        initializeTensorOperator();
    }
    virtual ~H_BM_dia() { }

protected:
    PositionOperator o_x;
    PositionOperator o_y;
    PositionOperator o_z;
    PositionOperator k_x;
    PositionOperator k_y;
    PositionOperator k_z;
    CubicPotential r_m3;

    void initializeTensorOperator() {
        static double alpha = 7.2973525664;

        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] = -(alpha*alpha/2.0)*r_m3*(o_y*k_y + o_z*k_z);
        h[0][1] =  (alpha*alpha/2.0)*r_m3*(o_x*k_y);
        h[0][2] =  (alpha*alpha/2.0)*r_m3*(o_x*k_z);
        h[1][0] =  (alpha*alpha/2.0)*r_m3*(o_y*k_x);
        h[1][1] = -(alpha*alpha/2.0)*r_m3*(o_x*k_x + o_z*k_z);
        h[1][2] =  (alpha*alpha/2.0)*r_m3*(o_y*k_z);
        h[2][0] =  (alpha*alpha/2.0)*r_m3*(o_z*k_x);
        h[2][1] =  (alpha*alpha/2.0)*r_m3*(o_z*k_y);
        h[2][2] = -(alpha*alpha/2.0)*r_m3*(o_x*k_x + o_y*k_y);
    }
};

#endif // H_BM_DIA_H

