//
// Created by machiels on 31/08/2022.
//

#include "ElementIntegrals.h"

ElementIntegrals::ElementIntegrals(double Rmid, double elemWidth) {
    double ratio = Rmid / elemWidth;

    //returns the integral of 1/R * phi_i * phi_j, aka I_{1,-1}. The first subscript is for the format, the second for the power of R in the integrand.
    moment1m1(0, 0) = (-2 * elemWidth * (elemWidth + Rmid) +
                       (elemWidth + 2 * Rmid) * (elemWidth + 2 * Rmid) * arccoth(2 * ratio)) /
                      (2 * elemWidth * elemWidth);
    moment1m1(0, 1) = ratio + (0.5 - 2 * ratio * ratio) * arccoth(2 * ratio);
    moment1m1(1,0) = moment1m1(0,1);
    moment1m1(1,1) = (2 * elemWidth * (elemWidth - Rmid) +
                      (elemWidth - 2 * Rmid) * (elemWidth - 2 * Rmid) * arccoth(2 * ratio)) /
                     (2 * elemWidth * elemWidth);

    //returns the integral of phi_i * phi_j, aka I_{1,0}
    moment1p0 << 2, 1, 1, 2;
    moment1p0 *= elemWidth / 6;

    //returns the integral of R * phi_i * phi_j, aka I_{1,1}
    moment1p1 << 4 * Rmid - elemWidth, 2 * Rmid, 2 * Rmid, 4 * Rmid + elemWidth;
    moment1p1 *= elemWidth / 12;

    //returns the integral of phi_i' * phi_j, aka I_{2,0}
    moment2p0 << -0.5, -0.5, 0.5, 0.5;

    //returns the integral of R * phi_i' * phi_j, aka I_{2,1}
    moment2p1 << elemWidth - 6 * Rmid, -elemWidth - 6 * Rmid, -elemWidth + 6 * Rmid, elemWidth + 6 * Rmid;
    moment2p1 /= 12;

    //returns the integral of phi_i * phi_j', aka I_{3,0}
    moment3p0 = moment2p0.transpose();

    //returns the integral of R * phi_i * phi_j', aka I_{3,1}
    moment3p1 = moment2p1.transpose();

    //returns the integral of R phi_i' * phi_j', aka I_{4,1}
    moment4p1 << 1, -1, -1, 1;
    moment4p1 *= ratio;

}

double ElementIntegrals::arccoth(double x) {
    return 0.5 * std::log((x + 1) / (x - 1));
}
