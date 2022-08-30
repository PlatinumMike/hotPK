//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_PLASMA_H
#define HOTPLASMAKERNEL_PLASMA_H

#include "Mesh.h"
#include "Species.h"

class Plasma {
public:
    Plasma(Mesh &mesh, double nTor, double omega);
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

    void addSpecies(double mass, double charge, double fraction, double omega);

private:
    Mesh *m_mesh;
    int m_NSpecies;
    std::vector<Species> specList{};
    double m_nTor;
    double m_omega;



};


#endif //HOTPLASMAKERNEL_PLASMA_H
