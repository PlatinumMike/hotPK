//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MATRIX_H
#define HOTPLASMAKERNEL_MATRIX_H

#include <complex>
#include <Eigen/Dense>
#include "Mesh.h"


class Matrix {
public:
    Matrix(int gridRes, Mesh &mesh, int nTor, double antFreq);

    /**
     * Element extraction routine
     * @param rowIndex row index
     * @param colIndex column index
     * @return value of the matrix (plasma contribution only)
     */
    std::complex<double> getEntry(int rowIndex, int colIndex);

    std::complex<double> getRhs(int rowIndex);

    /**
     * fills matrix and rhs
     */
    void buildMatrix();

    void solve();

private:
    const int m_gridRes;
    const int m_NDOF;
    const double m_omega;
    const int m_nTor;

    Eigen::MatrixXcd *globalMatrix;
    Eigen::VectorXcd *rhs;
    Eigen::VectorXcd *solution;

    Mesh *m_mesh;

    /**
     * converts global index to index on the mesh.
     * Using block ordering, so sorting by component first, then by node index.
     * @param i matrix row (or column) index
     * @return corresponding node index
     */
    int global2Node(int i);
    /**
     * Converts global index to index of potential (so 0,1,2 or 3)
     * @param i matrix row (or column) index
     * @return corresponding component index
     */
    int global2Comp(int i);

    ~Matrix();

};


#endif //HOTPLASMAKERNEL_MATRIX_H
