//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MATRIX_H
#define HOTPLASMAKERNEL_MATRIX_H

#include <complex>
#include <Eigen/Dense>
#include "Mesh.h"
#include "Plasma.h"
#include <string>



class Matrix {
public:
    Matrix(int gridRes, Mesh &mesh, Plasma &plasma, int nTor, double omega, plasmaType pType);

    /**
     * Element extraction routine.
     * Note: this method can be called directly after calling the constructor, no need to build the full matrix first.
     * @param rowIndex row index
     * @param colIndex column index
     * @return value of the matrix (plasma contribution only)
     */
    std::complex<double> getEntry(int rowIndex, int colIndex);

    std::complex<double> getRhs(int rowIndex);

    /**
     * Fills matrix and rhs. This is done by looping over all row and column indices.
     * Of course this is rather wasteful with the vacuum, cold and warm models, as those are local.
     * So the matrix will be sparse, with O(N) non-zeros. This construction takes O(N^2) of work.
     * However, the end goal is constructing the matrix for a hot plasma, where the matrix will in general be dense.
     * So for simplicity the same procedure is used.
     */
    void buildMatrix();

    void solve();

    void saveMatrix(std::string fileName, bool isReal);


private:
    const int m_gridRes;
    const int m_NDOF;
    const double m_omega;
    const int m_nTor;
    const plasmaType m_pType;

    Eigen::MatrixXcd *globalMatrix;
    Eigen::VectorXcd *rhs;
    Eigen::VectorXcd *solution;

    Mesh *m_mesh;
    Plasma *m_plasma;

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

    std::complex<double> getEntryPlasma(int rowIndex, int colIndex);
    std::complex<double> getEntryVacuum(int rowIndex, int colIndex);

    ~Matrix();

};


#endif //HOTPLASMAKERNEL_MATRIX_H
