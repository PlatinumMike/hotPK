//
// Created by machiels on 29/09/2022.
//

#ifndef HOTPLASMAKERNEL_PARAMETERS_H
#define HOTPLASMAKERNEL_PARAMETERS_H

struct Parameters{
    //bounds of domain
    double R_west; //left wall
    double R_east; //right wall
    double R_center; //middle
    double R0; //major radius of the magnetic axis, need not be the same as R_center
    double R_ant; //antenna position

    //grid size
    int gridResolution; //number of nodes.
    int numDegFreedom; //number of degrees of freedom, for linear FEM this is just 4 times the grid resolution. 4 because there is one DOF per node for every component of the potential.

    //other
    int nToroidal; //toroidal mode number
    double B0; //magnetic field strength at position R0
    double omega; //angular frequency
    double antennaPeriod;
    double omegaOverC;
    double J0R; //antenna scaling parameter, in R, Phi and Z directions
    double J0Phi;
    double J0Z;
    bool saveMatrix; //set to true to write the dense matrix to a csv file.
    bool solveMatrix; //set to true to solve matrix problem.

    //plasma
    double electronDensity; //peak density
    double minElecDensity; //minimum density
    double hydrogenConcentration; //hydrogen fraction, so n_H/n_e.
    double electronTemp; //peak electron temperature (J).
    double minTemp; //minimum temperature > 0 to avoid div by 0.
    int method; //toggle between cold, approx and full.
    int HarmonicMin; //minimal cyclotron harmonic to include //todo: not yet used. Also add the peaking factor as an input
    int HarmonicMax; //maximal cyclotron harmonic to include
    int angularResolution; //resolution used in evaluating the integrals for the plasma response.

    //todo: add method to print all of the loaded values inside of the code, so a user can confirm they are loaded correctly.
};

#endif //HOTPLASMAKERNEL_PARAMETERS_H
