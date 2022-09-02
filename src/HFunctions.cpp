//
// Created by machiels on 31/08/2022.
//

#include "HFunctions.h"

inline int getIndex(int ix, int iy) {
    constexpr int numHarmonics = 4; //0,1,2,3 (and negative)
    return ix * numHarmonics + iy;
}


inline constexpr double pi = boost::math::constants::pi<double>();
inline constexpr double sqrtPi = boost::math::constants::root_pi<double>();
inline constexpr double sqrtPiInv = boost::math::constants::one_div_root_pi<double>();

double HFunctions::getH(double x, int index, int harmonic, bool jacobian) {

    int signFlip = 1;
    if (harmonic<0 && (index==3 || index==5)){
        signFlip = -1;
    }
    int nAbs = std::abs(harmonic); //H_-n = H_n, except for H_2, H_4!

    int a0local = a0H[index];
    int a1local = a1H[index];
    if(jacobian){
        a0local++;
        a1local++;
    }

    double poly1 = (b0H[getIndex(index,nAbs)] + b1H[getIndex(index,nAbs)]*x*x + b2H[getIndex(index,nAbs)]*powBySquaring(x,4))*sqrtPiInv;
    double poly2 = c0H[getIndex(index,nAbs)] + c1H[getIndex(index,nAbs)]*x*x + c2H[getIndex(index,nAbs)]*powBySquaring(x,4) + c3H[getIndex(index,nAbs)]*powBySquaring(x,6);

    return signFlip*(powBySquaring(x,a0local) * poly1 * std::exp(-x*x) + powBySquaring(x,a1local) * poly2 * std::erfc(x));
}

double HFunctions::getU(double x, int index, int harmonic) {

    int signFlip = 1;
    if (harmonic<0 && (index==3 || index==5)){
        signFlip = -1;
    }
    int nAbs = std::abs(harmonic);

    int a0local = a0U[index];
    int a1local = a1U[index];

    double poly1 = (b0U[getIndex(index,nAbs)] + b1U[getIndex(index,nAbs)]*x*x + b2U[getIndex(index,nAbs)]*powBySquaring(x,4) + b3U[getIndex(index,nAbs)]*powBySquaring(x,6))*sqrtPiInv;
    double poly2 = c0U[getIndex(index,nAbs)] + c1U[getIndex(index,nAbs)]*x*x + c2U[getIndex(index,nAbs)]*powBySquaring(x,4) + c3U[getIndex(index,nAbs)]*powBySquaring(x,6) + c4U[getIndex(index,nAbs)]*powBySquaring(x,8);

    return signFlip * (powBySquaring(x,a0local) * poly1 * std::exp(-x*x) + powBySquaring(x,a1local) * poly2 * std::erfc(x));
}

double HFunctions::getV(double x, int index, int harmonic) {
    //copy-paste from Mathematica's CForm[]. Not per se optimized for rapid evaluation, but that can be improved later on.

    int signFlip = 1;
    if (harmonic < 0 && (index == 3 || index == 5)) {
        signFlip = -1;
    }
    int nAbs = std::abs(harmonic);

    if (index == 0) {
        if (nAbs == 0) {
            return (1 + 2 * x * x) / (16. * sqrtPi) * std::exp(-x * x);
        } else if (nAbs == 1) {
            return (-1 + 2 * x * x) / (48. * sqrtPi) * std::exp(-x * x) - (powBySquaring(x, 3) * std::erfc(x)) / 6.;
        } else if (nAbs == 2) {
            return (-1 + 14 * x * x + 192 * powBySquaring(x, 4)) / (240. * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 3) * (5 + 12 * x * x) * std::erfc(x)) / 15.;
        } else if (nAbs == 3) {
            return (-1 + 34 * x * x + 1312 * powBySquaring(x, 4) + 960 * powBySquaring(x, 6)) / (560. * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 3) * (35 + 224 * x * x + 120 * powBySquaring(x, 4)) * std::erfc(x)) / 70.;
        }
    } else if (index == 1) {
        if (nAbs == 0) {
            return -0.0625 * (3 + 2 * x * x) / (sqrtPi) * std::exp(-x * x);
        } else if (nAbs == 1) {
            return (1 + 2 * x * x) / (16. * sqrtPi) * std::exp(-x * x);
        } else if (nAbs == 2) {
            return -0.0125 * (-1 - 6 * x * x + 32 * powBySquaring(x, 4)) / (sqrtPi) * std::exp(-x * x) +
                   (2 * powBySquaring(x, 5) * std::erfc(x)) / 5.;
        } else if (nAbs == 3) {
            return -0.0017857142857142857 * (-3 - 38 * x * x + 576 * powBySquaring(x, 4) + 640 * powBySquaring(x, 6)) / (sqrtPi) *
                   std::exp(-x * x) + (8 * powBySquaring(x, 5) * (7 + 5 * x * x) * std::erfc(x)) / 35.;
        }
    } else if (index == 2) {
        if (nAbs == 0) {
            return 0;
        } else if (nAbs == 1) {
            return -0.06666666666666667 * (1 + x * x + 3 * powBySquaring(x, 4)) / (sqrtPi) * std::exp(-x * x) +
                   (powBySquaring(x, 3) * (5 + 6 * x * x) * std::erfc(x)) / 30.;
        } else if (nAbs == 2) {
            return (-2 * (-2 - 2 * x * x + 69 * powBySquaring(x, 4) + 30 * powBySquaring(x, 6))) / (105. * sqrtPi) * std::exp(-x * x) +
                   (powBySquaring(x, 3) * (35 + 168 * x * x + 60 * powBySquaring(x, 4)) * std::erfc(x)) / 105.;
        } else if (nAbs == 3) {
            return -0.0031746031746031746 * (-3 - 3 * x * x + 1101 * powBySquaring(x, 4) + 1480 * powBySquaring(x, 6) + 280 * powBySquaring(x, 8)) /
                   (sqrtPi) * std::exp(-x * x) +
                   (powBySquaring(x, 3) * (315 + 3402 * x * x + 3240 * powBySquaring(x, 4) + 560 * powBySquaring(x, 6)) * std::erfc(x)) / 630.;
        }
    } else if (index == 3) {
        if (nAbs == 0) {
            return 0;
        } else if (nAbs == 1) {
            return signFlip * (-0.1 * ((-1 + x) * (1 + x) * (1 + 2 * x * x)) / (sqrtPi) * std::exp(-x * x) +
                               (powBySquaring(x, 5) * std::erfc(x)) / 5.);
        } else if (nAbs == 2) {
            return signFlip * (-0.02857142857142857 * (1 + x * x + 18 * powBySquaring(x, 4) + 20 * powBySquaring(x, 6)) / (sqrtPi) *
                               std::exp(-x * x) + (4 * powBySquaring(x, 5) * (7 + 5 * x * x) * std::erfc(x)) / 35.);
        } else if (nAbs == 3) {
            return signFlip *
                   (-0.0015873015873015873 * (3 + 3 * x * x + 474 * powBySquaring(x, 4) + 1880 * powBySquaring(x, 6) + 560 * powBySquaring(x, 8)) /
                    (sqrtPi) * std::exp(-x * x) +
                    (powBySquaring(x, 5) * (567 + 1080 * x * x + 280 * powBySquaring(x, 4)) * std::erfc(x)) / 315.);
        }
    } else if (index == 4) {
        if (nAbs == 0) {
            return -0.125 * 1 / (sqrtPi) * std::exp(-x * x);
        } else if (nAbs == 1) {
            return (1 + 4 * x * x) / (24. * sqrtPi) * std::exp(-x * x) - (powBySquaring(x, 3) * std::erfc(x)) / 6.;
        } else if (nAbs == 2) {
            return ((1 + 4 * x * x) * (1 + 12 * x * x)) / (120. * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 3) * (5 + 6 * x * x) * std::erfc(x)) / 15.;
        } else if (nAbs == 3) {
            return (1 + 36 * x * x + 368 * powBySquaring(x, 4) + 160 * powBySquaring(x, 6)) / (280. * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 3) * (35 + 112 * x * x + 40 * powBySquaring(x, 4)) * std::erfc(x)) / 70.;
        }
    } else if (index == 5) {
        if (nAbs == 0) {
            return 0;
        } else if (nAbs == 1) {
            return signFlip * ((x * (-1 + 2 * x * x)) / (8. * sqrt2 * sqrtPi) * std::exp(-x * x) -
                               ((1 + 4 * powBySquaring(x, 4)) * std::erfc(x)) / (16. * sqrt2));
        } else if (nAbs == 2) {
            return signFlip * ((sqrt2 / sqrtPi * powBySquaring(x, 3) * (1 + x * x)) / 3 * std::exp(-x * x) -
                               (powBySquaring(x, 4) * (3 + 2 * x * x) * std::erfc(x)) / (3. * sqrt2));
        } else if (nAbs == 3) {
            return signFlip *
                   ((powBySquaring(x, 3) * (2 + 7 * x * x + 2 * powBySquaring(x, 4))) / (2. * sqrt2 * sqrtPi) * std::exp(-x * x) -
                    (powBySquaring(x, 4) * (9 + 16 * x * x + 4 * powBySquaring(x, 4)) * std::erfc(x)) / (4. * sqrt2));
        }
    } else if (index == 6) {
        if (nAbs == 0) {
            return -0.25 * x / (sqrt2 * sqrtPi) * std::exp(-x * x) - std::erfc(x) / (8. * sqrt2);
        } else if (nAbs == 1) {
            return (x * (1 + 2 * x * x)) / (8. * sqrt2 * sqrtPi) * std::exp(-x * x) -
                   ((-1 + 2 * x * x) * (1 + 2 * x * x) * std::erfc(x)) / (16. * sqrt2);
        } else if (nAbs == 2) {
            return (powBySquaring(x, 3) * (1 + 4 * x * x)) / (6. * sqrt2 * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 4) * (3 + 4 * x * x) * std::erfc(x)) / (6. * sqrt2);
        } else if (nAbs == 3) {
            return (powBySquaring(x, 3) * (1 + 13 * x * x + 6 * powBySquaring(x, 4))) / (6. * sqrt2 * sqrtPi) * std::exp(-x * x) -
                   (powBySquaring(x, 4) * (9 + 32 * x * x + 12 * powBySquaring(x, 4)) * std::erfc(x)) / (12. * sqrt2);
        }
    } else {
        return -1; //error: index not recognized
    }
    return -2; //error: harmonic not recognized
}

double HFunctions::powBySquaring(double x, int n) {
    //https://en.wikipedia.org/wiki/Exponentiation_by_squaring
    if (n < 0) {
        return powBySquaring(1 / x, -n);
    } else if (n == 0) {
        return 1;
    } else if (n % 2 == 0) {
        return powBySquaring(x * x, n / 2);
    } else {
        return x * powBySquaring(x * x, (n - 1) / 2);
    }
}
