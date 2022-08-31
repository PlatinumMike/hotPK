//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MATRIX_H
#define HOTPLASMAKERNEL_MATRIX_H

#include <complex>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "Mesh.h"
#include "Plasma.h"
#include <string>

typedef Eigen::Triplet<std::complex<double>> T;

class Matrix {
public:
    Matrix(int gridRes, Mesh &mesh, Plasma &plasma, int nTor, double omega);

    /**
     * Element extraction routine.
     * @param rowIndex row index
     * @param colIndex column index
     * @return value of the matrix (plasma contribution only)
     * @note this method can be called directly after calling the constructor, no need to build the full matrix first.
     * But the way the code is currently set-up, the sparse contribution of the matrix must be built first, but this is
     * O(N) amount of work, so is very fast.
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

    void saveMatrix(const std::string& fileName, bool isReal);

    void saveSolution(const std::string& fileName);


private:
    const int m_gridRes;
    const int m_NDOF;
    const double m_omega;
    const double m_omegaOverC2;
    const int m_nTor;

    Eigen::MatrixXcd *globalMatrix;
    Eigen::SparseMatrix<std::complex<double>> *sparseMatrix;
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
    int global2Node(int i) const;
    /**
     * Converts global index to index of potential (so 0,1,2 or 3)
     * @param i matrix row (or column) index
     * @return corresponding component index
     */
    int global2Comp(int i) const;

    /**
     * Get global index from the indices local to the element
     * @param elemIndex element index
     * @param localNodeIndex local index (0 or 1)
     * @param component quantity: Pot,AR,Aphi,AZ, so (0,1,2,3)
     * @return global index
     */
    int local2Global(int elemIndex, int localNodeIndex, int component);

    std::complex<double> getEntryPlasma(int rowIndex, int colIndex);

    /**
     * Retrieve value from the vacuum contribution (Laplace type operator) to the matrix
     * @param rowIndex  row index of the matrix
     * @param colIndex row index of the matrix
     * @return matrix entry
     * @warning The sparse matrix must be built before you use this element extraction routine!
     */
    std::complex<double> getEntryVacuum(int rowIndex, int colIndex);

    void buildRHS();
    void buildSparse();

    /**
     * Insert mini matrices into the larger local matrix.
     * @param localMatrix 8x8 block matrix
     * @param miniMatrix 2x2 block to be inserted
     * @param comp1 component1 corresponds to the test function G, FR, Fphi, FZ
     * @param comp2 component2 corresponds to the component of the solution: Pot, AR, Aphi or AZ.
     */
    void commitMiniMatrix(Eigen::Matrix<std::complex<double>, 8, 8> &localMatrix, const Eigen::Matrix2cd &miniMatrix,
                          int comp1, int comp2);

    /**
     * Add local matrix into the global matrix, well, actually into a triplet list which will sum into the global matrix.
     * @param tripletList
     * @param localMatrix 8x8 block to be inserted
     * @param elemIndex index of the current element
     */
    void commitLocalMatrix(std::vector<T> & tripletList, const Eigen::Matrix<std::complex<double>, 8, 8> &localMatrix, int elemIndex);

    ~Matrix();

};


#endif //HOTPLASMAKERNEL_MATRIX_H
