#ifndef PMRF_PHYSICSCONSTANTS_H
#define PMRF_PHYSICSCONSTANTS_H

namespace physConstants{
    //constants of nature
    constexpr double elementaryCharge = 1.602176634e-19; //elmentary charge (C), exact
    constexpr double massElectron = 9.1093837015e-31; //electron mass (kg)
    constexpr double massProton = 1.67262192369e-27; //proton mass (kg)
    constexpr double speedOfLight = 299792458.0; //speed of light (m/s), exact
    constexpr double mu_0 = 1.25663706212e-6; //vacuum magnetic permeability (H/m)
    constexpr double eps_0 = 8.8541878128e-12; //vacuum electric permittivity (F/m)
    constexpr double coulombConstant = 8.9875517923e9; //1/(4 pi eps_0) in kg m^3 s^-2 C^-2
    constexpr double mu_0inv = 1.0 / mu_0;
}

#endif //PMRF_PHYSICSCONSTANTS_H
