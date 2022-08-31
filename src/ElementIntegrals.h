//
// Created by machiels on 31/08/2022.
//

#ifndef HOTPLASMAKERNEL_ELEMENTINTEGRALS_H
#define HOTPLASMAKERNEL_ELEMENTINTEGRALS_H

#include <Eigen/Dense>

class ElementIntegrals {
public:
    ElementIntegrals(double Rmid, double elemWidth);
    //assemble local matrices, subscripts 'm' for minus, 'p' for plus.
    //todo: There is some precision loss due to cancellation issues (related to moment1m1 only), so recast this moment into a different form.
    Eigen::Matrix2d moment1m1, moment1p0, moment1p1, moment2p0, moment2p1, moment3p0, moment3p1, moment4p1;
private:

    /**
     * Computes inverse hyperbolic cotangent
     * @param x
     * @return arccoth(x)
     */
    static double arccoth(double x);

};


#endif //HOTPLASMAKERNEL_ELEMENTINTEGRALS_H
