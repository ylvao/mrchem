#ifndef H_BB_DIA_H
#define H_BB_DIA_H

#include "QMTensorOperator.h"
#include "H_E_dip.h"

class H_BB_dia : public RankTwoTensorOperator<3,3> {
public:
    H_BB_dia(const double *o = 0)
            : r_x(0, o), r_y(1, o), r_z(2, o) {
        initializeTensorOperator();
    }
    virtual ~H_BB_dia() { }

protected:
    PositionOperator r_x;
    PositionOperator r_y;
    PositionOperator r_z;

    void initializeTensorOperator() {
        RankTwoTensorOperator<3,3> &h = (*this);
        h[0][0] =  1.0/4.0*(r_y*r_y + r_z*r_z);
        h[0][1] = -1.0/4.0*(r_x*r_y);
        h[0][2] = -1.0/4.0*(r_x*r_z);
        h[1][0] = -1.0/4.0*(r_y*r_x);
        h[1][1] =  1.0/4.0*(r_x*r_x + r_z*r_z);
        h[1][2] = -1.0/4.0*(r_y*r_z);
        h[2][0] = -1.0/4.0*(r_z*r_x);
        h[2][1] = -1.0/4.0*(r_z*r_y);
        h[2][2] =  1.0/4.0*(r_x*r_x + r_y*r_y);
    }
};

#endif // H_BB_DIA_H

