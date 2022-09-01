//
// Created by machiels on 30/08/2022.
//

#include "Species.h"
#include "physicsConstants.h"
#include "gsl/gsl_sf_dawson.h"
#include "HFunctions.h"

constexpr double pi = 3.141592653589793;

Species::Species(double mass, double charge, double fraction, int nTor, double omegaIn, double temp, plasmaType pType) :
        m_mass(mass), m_charge(charge), m_fraction(fraction), m_nTor(nTor), omega(omegaIn), m_peakTemp(temp), m_pType(pType) {
    tempOffset = 1.0e2*physConstants::elementaryCharge; //minimal temperature

    //setup maps
    //HDict takes in a label, and returns the required H function for that kernel element
    HDict["R0"] = 0;
    HDict["R0mod"] = 1;
    HDict["R1"] = 2;
    HDict["R2"] = 3;
    HDict["R3"] = 4;
    HDict["R4"] = 5;
    HDict["R5"] = 6;
    HDict["D3"] = 4;
    HDict["D4"] = 5;
    HDict["D5"] = 6;

    //SDict takes in a label, and returns the required S function for that kernel element
    SDict["R0"] = 0;
    SDict["R0mod"] = 0;
    SDict["R1"] = 0;
    SDict["R2"] = 0;
    SDict["R3"] = 2;
    SDict["R4"] = 1;
    SDict["R5"] = 1;
    SDict["D3"] = 1;
    SDict["D4"] = 0;
    SDict["D5"] = 0;

    //some left over coefficients
    coefDict["R0"] = 2.0;
    coefDict["R0mod"] = -2.0;
    coefDict["R1"] = 1.0;
    coefDict["R2"] = {0, -1};
    coefDict["R3"] = -1.0;
    coefDict["R4"] = {0, -0.5 * std::sqrt(2) * m_charge / std::abs(m_charge)};
    coefDict["R5"] = -0.5 * std::sqrt(2) * m_charge / std::abs(m_charge);

}
std::complex<double> Species::getConductivity(double R, int row, int col) const {
    if(m_pType==vacuum){
        return {0};
    }
    //else: see below

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
    return -2.0 * gsl_sf_dawson(x) + std::sqrt(pi) * std::complex<double>{0.0, 1.0} * std::exp(-x * x);
}

std::complex<double> Species::plasmaDispersionFunctionDeriv(double x) const{
    return -2.0 * (1.0 + x * plasmaDispersionFunction(x));
}



std::complex<double> Species::getKernel(double R, int row, int col, double s1, double s2, double angle, int power) const {
    double cyc = getCyclotronFreq(R);
    double vT = getVThermal(R);
    std::complex<double> ans{0};
    double kpar = m_nTor/R;

    //todo: any way to do this more efficiently/cleanly? Not all elements are needed.
    std::complex<double> D3 = getElement("D3", R, s1, s2, kpar, power);
    std::complex<double> D4 = getElement("D4", R, s1, s2, kpar, power);
    std::complex<double> D5 = getElement("D5", R, s1, s2, kpar, power);
    std::complex<double> R0 = getElement("R0", R, s1, s2, kpar, power);
    std::complex<double> R0mod = getElement("R0mod", R, s1, s2, kpar, power);
    std::complex<double> R1 = getElement("R1", R, s1, s2, kpar, power);
    std::complex<double> R2 = getElement("R2", R, s1, s2, kpar, power);
    std::complex<double> R3 = getElement("R3", R, s1, s2, kpar, power);
    std::complex<double> R4 = getElement("R4", R, s1, s2, kpar, power);
    std::complex<double> R5 = getElement("R5", R, s1, s2, kpar, power);

    if (row == 0) {
        if (col == 0) {
            ans = -D4 * std::cos(angle) - D5 * std::sin(angle);
        } else if (col == 1) {
            ans = R1 + 0.5 * (R0 - R0mod * std::cos(2 * angle));
        } else if (col == 2) {
            ans = R2 - 0.5 * R0mod * std::sin(2 * angle);
        } else {
            ans = R4 * std::cos(angle) + R5 * std::sin(angle);
        }
    } else if (row == 1) {
        if (col == 0) {
            ans = -D4 * std::sin(angle) + D5 * std::cos(angle);
        } else if (col == 1) {
            ans = -R2 - 0.5 * R0mod * std::sin(2 * angle);
        } else if (col == 2) {
            ans = R1 + 0.5 * (R0 + R0mod * std::cos(2 * angle));
        } else {
            ans = R4 * std::sin(angle) - R5 * std::cos(angle);
        }
    } else {
        if (col == 0) {
            ans = -D3;
        } else if (col == 1) {
            ans = R4 * std::cos(angle) - R5 * std::sin(angle);
        } else if (col == 2) {
            ans = R4 * std::sin(angle) + R5 * std::cos(angle);
        } else {
            ans = R3;
        }
    }


    return 4.0 * pi * physConstants::eps_0 * omega * getPlasmaFreq2(R) * cyc * cyc / (vT * vT * vT) * ans;
}

double Species::thermalLarmorRadius(double R) const {
    return getVThermal(R)/std::abs(getCyclotronFreq(R));
}

double Species::getPerpIntegral(double xi1, double xi2, int indexH, int harmonic, double rhoT, int power) const{
    if (power == 1) {
        //the 4\rho^2 comes from the change of variables from s to \xi in the integration.
        return 4 * rhoT * rhoT * (HFunctions::getU(xi2, indexH, harmonic) - HFunctions::getU(xi1, indexH, harmonic));
    } else if (power == 2) {
        return 8 * rhoT * rhoT * rhoT * (HFunctions::getV(xi2, indexH, harmonic) - HFunctions::getV(xi1, indexH, harmonic));
    } else {
        return -1; //todo: other powers are not yet supported, throw error.
    }
}

std::complex<double> Species::getSIntegral(double R, int index, int harmonic, double kpar) const{
    double zeta = getZeta(R, harmonic);
    if (index == 0) {
        return 2 * pi / std::abs(kpar) * plasmaDispersionFunction(zeta);
    } else if (index == 1) {
        return 2 * pi / kpar * plasmaDispersionFunctionDeriv(zeta);
    } else {
        return 2 * pi / std::abs(kpar) * zeta * plasmaDispersionFunctionDeriv(zeta);
    }
}

std::complex<double> Species::getElement(const std::string& label, double R, double s1, double s2, double kpar, int power) const{
    int indexH = HDict.at(label);
    int indexS = SDict.at(label);
    double rhoT = thermalLarmorRadius(R);
    double xi1 = s1 / (2 * rhoT);
    double xi2 = s2 / (2 * rhoT);
    std::complex<double> ans{0};
    //for now loop over harmonics -3 to +3, which should be enough for nearly all scenarios, but consider extending it in the future.
    for (int harmonic = -3; harmonic < 4; harmonic++) {
        ans += getPerpIntegral(xi1, xi2, indexH, harmonic, rhoT, power) * getSIntegral(R, indexS, harmonic, kpar);
    }
    return getCoef(label, R) * ans;
}

std::complex<double> Species::getCoef(const std::string& label, double R) const {
    if (label == "D3") {
        return -physConstants::speedOfLight / getVThermal(R);
    } else if (label == "D4") {
        return {0, std::sqrt(2) * m_charge / std::abs(m_charge) * physConstants::speedOfLight / getVThermal(R)};
    } else if (label == "D5") {
        return std::sqrt(2) * m_charge / std::abs(m_charge) * physConstants::speedOfLight / getVThermal(R);
    } else {
        return coefDict.at(label);
    }
}
