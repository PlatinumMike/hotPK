//
// Created by machiels on 30/08/2022.
//

#ifndef HOTPLASMAKERNEL_SPECIES_H
#define HOTPLASMAKERNEL_SPECIES_H

#include <complex>
#include <map>

enum plasmaType {vacuum, cold, warm, hot};


class Species {
public:
    Species(double mass, double charge, double peakDensity, double minDens, int nTor, double omegaIn,
            double peakTemp, double minTemp, double RAxis, double RWest, double REast, double BAxis, plasmaType pType);

    /**
     * compute given element of the 3x3 conductivity kernel.
     * @param R major radius
     * @param row row index
     * @param col column index
     * @param species species index
     * @return \f$\sigma_{ij}\f$
     */
    [[nodiscard]] std::complex<double> getConductivity(double R, int row, int col) const;

    /**
     * Computes integral of kernel, times s^power. Already integrated over parallel direction
     * @param R major radius
     * @param row row index
     * @param col column index
     * @param s1 lower integration bound s
     * @param s2 upper integration bound s
     * @param angle angle
     * @param power power of s in integrand (1 or 2)
     * @return integral
     */
    std::complex<double> getKernel(double R, int row, int col, double s1, double s2, double angle, int power) const;

private:
    const double B0Axis;
    const double m_RAxis;
    const double m_REast;
    const double m_RWest;
    const double m_mass;
    const double m_charge;
    const double m_chargeSign;
    const double m_peakDensity;
    const double densOffset;
    const int m_nTor;
    const double omega;
    const plasmaType m_pType;
    const double m_peakTemp;
    const double tempOffset; //minimal temperature

    //helper maps for looping up indices
    std::map<std::string,int> HDict;
    std::map<std::string,int> SDict;
    std::map<std::string,std::complex<double>> coefDict;
    /**
     * Get coefficient needed for computing conductivity kernel
     * @param label element label
     * @param R major radius
     * @return coefficient
     */
    std::complex<double> getCoef(const std::string& label,double R) const;

    /**
     * Get S entry in the cold dielectric tensor
     * @param R major radius
     * @return \f$ S \f$
     */
    double getSCold(double R) const;
    /**
      * Get D entry in the cold dielectric tensor
      * @param R major radius
      * @return \f$ D \f$
      */
    double getDCold(double R) const;
    /**
      * Get P entry in the cold dielectric tensor
      * @param R major radius
      * @return \f$ P \f$
      */
    double getPCold(double R) const;

    // warm versions of the above
    std::complex<double> getSWarm(double R) const;
    std::complex<double> getDWarm(double R) const;
    std::complex<double> getPWarm(double R) const;

    /**
     * Get square of plasma frequency
     * @param R major radius
     * @return \f$ \omega_p^2 \f$
     */
    double getPlasmaFreq2(double R) const;
    /**
     * Get background magnetic field strength
     * @param R major radius
     * @return magnetic field strength
     */
    double getB0(double R) const;
    /**
     * Get sample from any profile (density/temperature)
     * @param R major radius
     * @param minValue minimal value, in the edge
     * @param maxValue peak value, in the core
     * @param peakingFactor decide how peaked the profile is, 1 means parabolic, higher is more peaked.
     * @return sample of profile
     */
    double getProfile(double R, double minValue, double maxValue, double peakingFactor) const;
    /**
     * Get species density
     * @param R major radius
     * @return density
     */
    double getDensity(double R) const;
    /**
     * Get cyclotron frequency
     * @param R major radius
     * @return frequency (rad/s)
     */
    double getCyclotronFreq(double R) const;
    /**
     * Get species temperature
     * @param R major radius
     * @return temperature (J)
     */
    double getTemperature(double R) const;
    /**
     * Get thermal velocity
     * @param R major radius
     * @return \f$ \sqrt{2 T/m}\f$ in (m/s)
     */
    double getVThermal(double R) const;
    /**
     * Get resonance parameter
     * @param R major radius
     * @param harmonic cyclotron harmonic. Integer value between -inf and +inf expected.
     * @return \f$ \zeta_n \f$
     */
    double getZeta(double R, int harmonic) const;

    /**
     * Compute plasma dispersion function using the Dawson function.
     * @param x input value
     */
    std::complex<double> plasmaDispersionFunction(double x) const;

    /**
     * Get derivative of plasma dispersion function
     * @param x input value
     * @return derivative of PDF
     */
    std::complex<double> plasmaDispersionFunctionDeriv(double x) const;

    /**
     * Compute Larmor radius using thermal velocity
     * @param R major radius
     * @return \f$ m v_T/(|q| B)  = v_T / |\Omega| \f$
     */
    double thermalLarmorRadius(double R) const;

    /**
     * Compute perpendicular integrals needed for hot plasma model
     * @param xi1 normalized lower bound radial distance
     * @param xi2 normalized upper bound radial distance
     * @param indexH index of the H function to be used
     * @param harmonic cyclotron harmonic
     * @param rhoT thermal Larmor radius
     * @param power power of s in the integrand
     * @return perpendicular integral
     */
    double getPerpIntegral(double xi1, double xi2, int indexH, int harmonic, double rhoT, int power) const;

    /**
     * Compute integral in parallel direction, this is done in the usual Fourier space actually, not configuration space.
     * @param R major radius
     * @param index index of S function to be used (well, it's Fourier transform)
     * @param harmonic cyclotron harmonic
     * @param kpar parallel wave number
     * @return parllel integral
     */
    std::complex<double> getSIntegral(double R, int index, int harmonic, double kpar) const;

    /**
     * Compute entry to conductivity kernel, integrated over the radial distance.
     * @param label name of the entry
     * @param R major radius
     * @param s1 lower integration bound radial distance
     * @param s2 upper integration bound radial distance
     * @param kpar parallel wave number
     * @param power power of s in the integrand
     * @return entry in conductivity kernel, integrated
     */
    std::complex<double> getElement(const std::string& label, double R, double s1, double s2, double kpar, int power) const;

};


#endif //HOTPLASMAKERNEL_SPECIES_H
