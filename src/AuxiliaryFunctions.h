//
// Created by machiels on 01/09/2022.
//

#ifndef HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H
#define HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H

#include <complex>

class AuxiliaryFunctions {
public:
    static void linspace(double *array, double start, double stop, int count);

    static void getAngle(double *angleArray, int resolution, int sector);

    /**
     * composit Simpson rule to integrate \int dx y(x)
     * @param x position array, size N
     * @param y values y(x)
     * @return integral
     * @warning N should be odd
     */
    template <typename T>
    static T integrateSimpson(const double *x, const T *y, int size) {
        //Composite Simpson's rule for irregularly spaced data
        //https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data
        T answer = 0;

        int intervals = size - 1; //intervals assumed to be even! //todo: handle odd case.
        double h2i, h2iplus1; //grid spacings, h_i = x_{i+1} - x_i
        for (int i = 0; i < intervals / 2; i++) {
            h2i = x[2 * i + 1] - x[2 * i];
            h2iplus1 = x[2 * i + 2] - x[2 * i + 1];
            answer += (h2i + h2iplus1) / 6 * ((2 - h2iplus1 / h2i) * y[2 * i]
                                              + (h2i + h2iplus1) * (h2i + h2iplus1) / (h2i * h2iplus1) * y[2 * i + 1]
                                              + (2 - h2i / h2iplus1) * y[2 * i + 2]);
        }
        return answer;
    }
};


#endif //HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H
