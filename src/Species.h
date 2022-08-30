//
// Created by machiels on 30/08/2022.
//

#ifndef HOTPLASMAKERNEL_SPECIES_H
#define HOTPLASMAKERNEL_SPECIES_H

#include <complex>
enum plasmaType {vacuum, cold, warm, hot};


class Species {
public:
    Species(double mass, double charge, double fraction, double omegaIn);

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
    const double B0Axis = 3; //todo: magic constants, fix
    const double RAxis = 3;
    const double m_mass;
    const double m_charge;
    const double m_fraction;
    const double peakElectronDensity = 1.0e20;
    const double omega;

    //dielectric tensor elements
    double getSCold(double R) const;
    double getDCold(double R) const;
    double getPCold(double R) const;

    double getPlasmaFreq2(double R) const;
    double getB0(double R) const;
    double getDensity(double R) const;
    double getCyclotronFreq(double R) const;

};


#endif //HOTPLASMAKERNEL_SPECIES_H
