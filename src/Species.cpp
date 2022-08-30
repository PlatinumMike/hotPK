//
// Created by machiels on 30/08/2022.
//

#include "Species.h"
#include "physicsConstants.h"

Species::Species(double mass, double charge, double fraction, double omegaIn) : m_mass(mass), m_charge(charge), m_fraction(fraction), omega(omegaIn){

}
std::complex<double> Species::getConductivity(double R, int row, int col) const {
    std::complex<double> ans{0};
    //the matrix ordering is determined from this: R->y, phi->z, Z->x
    //for the cold plasma some entries are 0, so only overwriting the non-zero ones.
    if (row == 0) {
        if (col == 0) {
            ans = getSCold(R);
        } else if (col == 2) {
            ans = std::complex<double>{0, 1} * getDCold(R);
        }
    } else if (row == 1) {
        if (col == 1) {
            ans = getPCold(R);
        }
    } else {
        if (col == 0) {
            ans = -std::complex<double>{0, 1} * getDCold(R);
        } else if (col == 2) {
            ans = getSCold(R);
        }
    }
    return -std::complex<double>{0, 1} * omega * physConstants::eps_0 * ans;
}

double Species::getSCold(double R) const{
    double cyc = getCyclotronFreq(R);
    //todo: add tiny im part in omega to avoid div by 0.
    return -getPlasmaFreq2(R) / (omega * omega - cyc * cyc);
}

double Species::getDCold(double R) const{
    double cyc = getCyclotronFreq(R);
    return cyc * getPlasmaFreq2(R) / (omega * (omega * omega - cyc * cyc));
}

double Species::getPCold(double R) const{
    return -getPlasmaFreq2(R)/(omega*omega);
}

double Species::getPlasmaFreq2(double R) const{
    return m_charge*m_charge/(physConstants::eps_0*m_mass)*getDensity(R);
}

double Species::getB0(double R) const{
    return B0Axis * RAxis / R;
}

double Species::getDensity(double R) const{
    return std::max(m_fraction * peakElectronDensity * (1.0 - 3.0 * (R - RAxis) * (R - RAxis)), 0.0);
}

double Species::getCyclotronFreq(double R) const{
    return m_charge* getB0(R)/m_mass;
}
