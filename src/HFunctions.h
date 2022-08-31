//
// Created by machiels on 31/08/2022.
//

#ifndef HOTPLASMAKERNEL_HFUNCTIONS_H
#define HOTPLASMAKERNEL_HFUNCTIONS_H

#include <array>
#include <cmath>
#include "boost/math/constants/constants.hpp"

class HFunctions {
public:
    /**
     * Compute H function
     * @param x position (normalized)
     * @param index select H0, H0breve, H1, H2, H3, H4, H5 with index 0,1,2,3,4,5,6
     * @param harmonic harmonic
     * @param jacobian optional, if true, multiply by x to avoid 1/x singularity. False by default.
     * @return H(x)
     * @warning harmonic only mapped in the range -3 to 3, no check is performed!
     */
    //todo: implement check on the index and harmonic so those do not go out of range!
    static double getH(double x, int index, int harmonic, bool jacobian=false);

    /**
     * Integral of H, \f$ \int dx x H(x) \f$
     * @param x x
     * @param index select U0, U0breve, U1, U2, U3, U4, U5
     * @param harmonic harmonic
     * @return \f$ U(x) \f$
     */
    static double getU(double x, int index, int harmonic);

    /**
    * Integral of H, \f$ \int dx x^2 H(x) \f$
    * @param x x
    * @param index select V0, V0breve, V1, V2, V3, V4, V5
    * @param harmonic harmonic
    * @return \f$ V(x) \f$
    */
    static double getV(double x, int index, int harmonic);

private:

    inline static constexpr double sqrt2 = boost::math::constants::root_two<double>();

    //constants needed for H(x). Note, this is not very scalable if you need to go to arbitrary harmonic, for that you are probably better off using the hypergeometric functions.
    static constexpr std::array<int,7> a0H = {-1,-1,1 ,1 ,-1,0 ,0};
    static constexpr std::array<int,7> a1H = {0 , 0, 0,0 ,0 ,1 ,1};

    static constexpr std::array<double, 28> b0H = {
            1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
            1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
            0, -1, -6, -15,
            0, -1, -2, -3,
            1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0,
            0, 1 / sqrt2, sqrt2, 3 / sqrt2,
            1 / (2 * sqrt2), 1 / (2 * sqrt2), 1 / (2 * sqrt2), 1 / (2 * sqrt2)};

    static constexpr std::array<double, 28> b1H = {-1.0 / 4.0, 1.0 / 4.0, 15.0 / 4.0, 41.0 / 4.0,
                                  1.0 / 4.0, -1.0 / 4.0, -7.0 / 4.0, -17.0 / 4.0,
                                  0, 0, -4, -32,
                                  0, 0, -4, -20,
                                  0, 0, 2, 6,
                                  0, 0, 2 * sqrt2, 10 * sqrt2,
                                  0, 0, 2 * sqrt2, 6 * sqrt2};

    static constexpr std::array<double, 28> b2H = {0, 0, 0, 12,
                                  0, 0, 0, -8,
                                  0, 0, 0, -8,
                                  0, 0, 0, -8,
                                  0, 0, 0, 4,
                                  0, 0, 0, 4 * sqrt2,
                                  0, 0, 0, 4 * sqrt2};

    static constexpr std::array<double, 28> c0H = {0, -0.5, -1, -1.5,
                                  0, 0, 0, 0,
                                  0, 0.5, 1, 1.5,
                                  0, 0, 0, 0,
                                  0, -0.5, -1, -1.5,
                                  0, -1 / sqrt2, -2 * sqrt2, -9 / sqrt2,
                                  0, -1 / sqrt2, -sqrt2, -3 / sqrt2};
    static constexpr std::array<double, 28> c1H = {0, 0, -4, -16, 0, 0, 2, 8, 0, 1, 8, 27, 0, 1, 4, 9,
                                  0, 0, -2, -8, 0, 0, -2 * sqrt2, -12 * sqrt2, 0, 0, -2 * sqrt2,
                                  -8 * sqrt2};
    static constexpr std::array<double, 28> c2H = {0, 0, 0, -12, 0, 0, 0, 8, 0, 0, 4, 36, 0, 0, 4, 24, 0, 0, 0, -4,
                                  0, 0, 0, -4 * sqrt2, 0, 0, 0, -4 * sqrt2};
    static constexpr std::array<double, 28> c3H = {0, 0, 0, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0, 8,
                                  0, 0, 0, 8,
                                  0, 0, 0, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0, 0};

    //constants needed for U(x)
    static constexpr std::array<int, 7> a0U = {1, 1, 1, 1, 1, 0, 0};
    static constexpr std::array<int, 7> a1U = {2, 0, 0, 0, 0, 3, 3};

    static constexpr std::array<double, 28> b0U = {0.125, 0.125, 0.125, 0.125,
                                                   -0.125, 0.125, 0.125, 0.125,
                                                   0, -0.125, 0, 0,
                                                   0, 0.125, 0, 0,
                                                   0, 0.25, 0.25, 0.25,
                                                   0, -1 / (6 * sqrt2), -1 / (15 * sqrt2), -3 / (70 * sqrt2),
                                                   -1 / (4 * sqrt2), 1 / (12 * sqrt2), 1 / (60 * sqrt2),
                                                   1 / (140 * sqrt2)};

    static constexpr std::array<double, 28> b1U = {0, 0, 1, 3,
                                                   0, 0, -0.5, -4.0 / 3.0,
                                                   0, -0.25, -5.0 / 3.0, -9.0 / 2.0,
                                                   0, -0.25, -2.0 / 3.0, -1,
                                                   0, 0, 0.5, 5.0 / 3.0,
                                                   0, 1 / (3 * sqrt2), 7 * sqrt2 / 15, 51 / (35 * sqrt2),
                                                   0, 1 / (3 * sqrt2), 2 * sqrt2 / 15, 9 / (35 * sqrt2)};

    static constexpr std::array<double, 28> b2U = {0, 0, 0, 2,
                                                   0, 0, 0, -4.0 / 3.0,
                                                   0, 0, -2.0 / 3.0, -11.0 / 2.0,
                                                   0, 0, -2.0 / 3.0, -7.0 / 2.0,
                                                   0, 0, 0, 2.0 / 3.0,
                                                   0, 0, 2 * sqrt2 / 5, 74 * sqrt2 / 35,
                                                   0, 0, 2 * sqrt2 / 5, 46 * sqrt2 / 35};

    static constexpr std::array<double, 28> b3U = {0, 0, 0, 0,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, -1,
                                                   0, 0, 0, -1,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, 4 * sqrt2 / 7,
                                                   0, 0, 0, 4 * sqrt2 / 7};

    static constexpr std::array<double, 28> c0U = {0, -0.25, -0.5, -3.0 / 4.0,
                                                   -0.125, 0, 0, 0,
                                                   0, -1.0 / 16.0, 0, 0,
                                                   0, 1.0 / 16.0, 0, 0,
                                                   -0.125, 0, 0, 0,
                                                   0, -1 / (3 * sqrt2), -2 * sqrt2 / 3, -3 / sqrt2,
                                                   0, -1 / (3 * sqrt2), -sqrt2 / 3, -1 / sqrt2};

    static constexpr std::array<double, 28> c1U = {0, 0, -1, -4,
                                                   0, 0, 0, 0,
                                                   0, 0.25, 0.5, 3.0 / 4.0,
                                                   0, 0, 0, 0,
                                                   0, -0.25, -0.5, -3.0 / 4.0,
                                                   0, 0, -2 * sqrt2 / 5, -12 * sqrt2 / 5,
                                                   0, 0, -2 * sqrt2 / 5, -8 * sqrt2 / 5};

    static constexpr std::array<double, 28> c2U = {0, 0, 0, -2,
                                                   0, 0, 0.5, 2,
                                                   0, 0.25, 2, 27.0 / 4.0,
                                                   0, 0.25, 1, 9.0 / 4.0,
                                                   0, 0, -0.5, -2,
                                                   0, 0, 0, -4 * sqrt2 / 7,
                                                   0, 0, 0, -4 * sqrt2 / 7};

    static constexpr std::array<double, 28> c3U = {0, 0, 0, 0,
                                                   0, 0, 0, 4.0 / 3.0,
                                                   0, 0, 2.0 / 3.0, 6,
                                                   0, 0, 2.0 / 3.0, 4,
                                                   0, 0, 0, -2.0 / 3.0,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, 0};

    static constexpr std::array<double, 28> c4U = {0, 0, 0, 0,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, 1,
                                                   0, 0, 0, 1,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, 0,
                                                   0, 0, 0, 0};

};


#endif //HOTPLASMAKERNEL_HFUNCTIONS_H
