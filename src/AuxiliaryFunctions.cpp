//
// Created by machiels on 01/09/2022.
//

#include "AuxiliaryFunctions.h"



void AuxiliaryFunctions::linspace(double *array, const double start, const double stop, const int count){
    double delta = (stop - start) / (count - 1);
    for (int i = 0; i < count; ++i) {
        array[i] = start + delta * i;
    }
}

void AuxiliaryFunctions::getAngle(double *angleArray, int resolution, int sector) {
    constexpr double pi = 3.141592653589793;
    if (sector == 1) {
        linspace(angleArray, 0.5 * pi, 1.5 * pi, resolution);
    } else if (sector == 2) {
        linspace(angleArray, -0.5 * pi, 1.5 * pi, resolution);
    } else {
        linspace(angleArray, -0.5 * pi, 0.5 * pi, resolution);
    }
}
