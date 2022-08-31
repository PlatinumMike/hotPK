//
// Created by machiels on 30/08/2022.
//

#ifndef HOTPLASMAKERNEL_SPECIES_H
#define HOTPLASMAKERNEL_SPECIES_H

#include <complex>
enum plasmaType {vacuum, cold, warm, hot};


class Species {
public:
    Species(double mass, double charge, double fraction, int nTor, double omegaIn, double temp, plasmaType pType);

    /**
     * compute given element of the 3x3 conductivity kernel.
     * @param R major radius
     * @param row row index
     * @param col column index
     * @param species species index
     * @return \f$\sigma_{ij}\f$
     */
    [[nodiscard]] std::complex<double> getConductivity(double R, int row, int col) const;

private:
    const double B0Axis = 2.623778994743588; //todo: magic constants, fix
    const double RAxis = 3;
    const double m_mass;
    const double m_charge;
    const double m_fraction;
    const double peakElectronDensity = 1.0e20;
    const int m_nTor;
    const double omega;
    const plasmaType m_pType;
    const double m_peakTemp;
    double tempOffset;

    //dielectric tensor elements
    double getSCold(double R) const;
    double getDCold(double R) const;
    double getPCold(double R) const;

    std::complex<double> getSWarm(double R) const;
    std::complex<double> getDWarm(double R) const;
    std::complex<double> getPWarm(double R) const;

    double getPlasmaFreq2(double R) const;
    double getB0(double R) const;
    double getDensity(double R) const;
    double getCyclotronFreq(double R) const;
    double getTemperature(double R) const;
    double getVThermal(double R) const;
    double getZeta(double R, int harmonic) const;

    /**
     * Compute plasma dispersion function using the Dawson function.
     * @param x input value
     */
    std::complex<double> plasmaDispersionFunction(double x) const;

    std::complex<double> plasmaDispersionFunctionDeriv(double x) const;

};


#endif //HOTPLASMAKERNEL_SPECIES_H
