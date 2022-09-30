//
// Created by machiels on 01/09/2022.
//

#ifndef HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H
#define HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H

#include <complex>

class AuxiliaryFunctions {
public:
    /**
     * Fill array with linearly spaced values
     * @param array target array
     * @param start starting value
     * @param stop final value
     * @param count number of entries to fill
     */
    static void linspace(double *array, double start, double stop, int count);

    /**
     * FIll array of polar angles needed for the integration
     * @param angleArray target array
     * @param resolution number of entries
     * @param sector select angular range, based on where we are w.r.t. the tent function (left, middle or right).
     */
    static void getAngle(double *angleArray, int resolution, int sector);

    /**
     * composit Simpson rule to integrate \int dx y(x)
     * @param x position array, size N
     * @param y values y(x)
     * @return integral
     */
    template <typename T>
    static T integrateSimpson(const double *x, const T *y, int size) {
        //Composite Simpson's rule for irregularly spaced data
        //https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data
        T answer = 0;

        int intervals = size - 1;
        double h2i, h2iplus1; //grid spacings, h_i = x_{i+1} - x_i
        for (int i = 0; i < intervals / 2; i++) {
            h2i = x[2 * i + 1] - x[2 * i];
            h2iplus1 = x[2 * i + 2] - x[2 * i + 1];
            answer += (h2i + h2iplus1) / 6 * ((2 - h2iplus1 / h2i) * y[2 * i]
                                              + (h2i + h2iplus1) * (h2i + h2iplus1) / (h2i * h2iplus1) * y[2 * i + 1]
                                              + (2 - h2i / h2iplus1) * y[2 * i + 2]);
        }

        if (intervals % 2 == 1) {
            //handle final interval separately
            h2i = x[intervals - 1] - x[intervals - 2]; //one before last interval
            h2iplus1 = x[intervals] - x[intervals - 1]; //last interval
            answer += h2iplus1 * (2 * h2iplus1 + 3 * h2i) / (6 * (h2i + h2iplus1)) * y[intervals];
            answer += h2iplus1 * (h2iplus1 + 3 * h2i) / (6 * h2i) * y[intervals - 1];
            answer -= h2iplus1 * h2iplus1 * h2iplus1 / (6 * h2i * (h2i + h2iplus1)) * y[intervals - 2];
        }

        return answer;
    }
};


#endif //HOTPLASMAKERNEL_AUXILIARYFUNCTIONS_H
