//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_PLASMA_H
#define HOTPLASMAKERNEL_PLASMA_H

#include "Mesh.h"
#include "Species.h"

class Plasma {
public:
    Plasma(Mesh &mesh, int nTor, double omega, int angularResolution, double REast, double RWest, plasmaType pType);
    /**
     * Gets 3x4 matrix that relates current to the 4 potentials.
     * Plasma response due to the basis function with index nodeIndex.
     * \f$J_i = \sum_{j,k} c_{ijk} A_{jk}\f$
     * @param R major radius
     * @param row row index i of the matrix, so the component of the current
     * @param col column index j of the matrix, so the component of the potential
     * @param nodeIndex index of the mesh node k
     * @return \f$c_{ijk}\f$
     */
    std::complex<double> getCurrentMatrix(double R, int row, int col, int nodeIndex);

    /**
     * Computes contribution of one species to current density matrix
     * The hot version is non-local, so requires species treatment.
     * @param R major radius
     * @param row row index
     * @param col column index
     * @param nodeIndex index of the mesh node
     * @param mesh mesh needed to evaluate the integrals
     * @return entry in 3x4 block matrix
     */
    std::complex<double> getCurrentMatrixHot(double R, int row, int col, int nodeIndex);

    void addSpecies(double mass, double charge, double peakDensity, double minDens, double omega, double peakTemp,
                    double minTemp, double R0, double B0, plasmaType pType);

private:
    Mesh *m_mesh;
    int m_NSpecies;
    std::vector<Species> specList{};
    int m_nTor;
    double m_omega;
    const int m_angularResolution;
    plasmaType m_pType;
    const double m_REast, m_RWest;



};


#endif //HOTPLASMAKERNEL_PLASMA_H
