//
// Created by machiels on 30/08/2022.
//

#include "Species.h"
#include "physicsConstants.h"
#include "gsl/gsl_sf_dawson.h"

Species::Species(double mass, double charge, double fraction, int nTor, double omegaIn, double temp, plasmaType pType) :
        m_mass(mass), m_charge(charge), m_fraction(fraction), m_nTor(nTor), omega(omegaIn), m_peakTemp(temp), m_pType(pType) {
    tempOffset = 1.0e2*physConstants::elementaryCharge; //minimal temperature

}
std::complex<double> Species::getConductivity(double R, int row, int col) const {
    if(m_pType==vacuum){
        return {0};
    }
    //else:

    std::complex<double> ans{0};
    //the matrix ordering is determined from this: R->y, phi->z, Z->x
    //for the cold, warm plasma some entries are 0, so only overwriting the non-zero ones.
    //todo: lot of code repetition here, can this be improved?
    if (m_pType == cold) {
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
    } else if (m_pType == warm) {
        if (row == 0) {
            if (col == 0) {
                ans = getSWarm(R);
            } else if (col == 2) {
                ans = std::complex<double>{0, 1} * getDWarm(R);
            }
        } else if (row == 1) {
            if (col == 1) {
                ans = getPWarm(R);
            }
        } else {
            if (col == 0) {
                ans = -std::complex<double>{0, 1} * getDWarm(R);
            } else if (col == 2) {
                ans = getSWarm(R);
            }
        }
    } else {
        return -1; //todo implement hot plasma
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

std::complex<double> Species::getSWarm(double R) const {
    double kpar = m_nTor / R; //only for a purely toroidal field.
    double zetaMinus1 = getZeta(R, -1);
    double zetaPlus1 = getZeta(R, 1);
    return getPlasmaFreq2(R) / (2 * omega * std::abs(kpar) * getVThermal(R)) *
           (plasmaDispersionFunction(zetaMinus1) + plasmaDispersionFunction(zetaPlus1));
}

std::complex<double> Species::getDWarm(double R) const {
    double kpar = m_nTor / R;
    double zetaMinus1 = getZeta(R, -1);
    double zetaPlus1 = getZeta(R, 1);
    return getPlasmaFreq2(R) / (2 * omega * std::abs(kpar) * getVThermal(R)) *
           (plasmaDispersionFunction(zetaMinus1) - plasmaDispersionFunction(zetaPlus1));
}

std::complex<double> Species::getPWarm(double R) const {
    double kpar = m_nTor / R;
    double zeta0 = getZeta(R, 0);
    return -getPlasmaFreq2(R)/(omega*std::abs(kpar)* getVThermal(R))*zeta0* plasmaDispersionFunctionDeriv(zeta0);
}

double Species::getTemperature(double R) const {
    return std::max(m_peakTemp * (1.0 - 3.0 * (R - RAxis) * (R - RAxis)), tempOffset);
}

double Species::getVThermal(double R) const {
    return std::sqrt(2 * getTemperature(R) / m_mass);
}

double Species::getZeta(double R, int harmonic) const {
    double kpar = m_nTor/R;
    return (omega - harmonic* getCyclotronFreq(R))/(std::abs(kpar)* getVThermal(R));
}

std::complex<double> Species::plasmaDispersionFunction(double x) const{
    //pdf = i \sqrt{\pi} erfc(ix), but complex arguments not supported currently by c++ standard library nor BOOST.
    //you can also write it as -2Dawson(x) + i \sqrt{\pi} \exp{-x^2}, to get the Dawson function the GNU Scientific Library is used.
    constexpr double pi = 3.141592653589793;
    return -2.0 * gsl_sf_dawson(x) + std::sqrt(pi) * std::complex<double>{0.0, 1.0} * std::exp(-x * x);
}

std::complex<double> Species::plasmaDispersionFunctionDeriv(double x) const{
    return -2.0 * (1.0 + x * plasmaDispersionFunction(x));
}
