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
    constexpr double forbiddenAngle = 1.0e-4; //just to avoid div by 0 when angle is near pi/2, 3pi/2, etc.
    if (sector == 1) {
        linspace(angleArray, 0.5 * pi + forbiddenAngle, 1.5 * pi - forbiddenAngle, resolution);
    } else if (sector == 2) {
        int halfRes = resolution / 2;
        linspace(angleArray, -0.5 * pi + forbiddenAngle, 0.5 * pi - forbiddenAngle, halfRes);
        linspace(angleArray + halfRes, 0.5 * pi + forbiddenAngle, 1.5 * pi - forbiddenAngle, resolution - halfRes);
    } else {
        linspace(angleArray, -0.5 * pi + forbiddenAngle, 0.5 * pi - forbiddenAngle, resolution);
    }
}
