//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MATRIX_H
#define HOTPLASMAKERNEL_MATRIX_H

#include <complex>


class Matrix {
public:
    Matrix(int gridRes);
    /**
     * Element extraction routine
     * @param rowIndex row index
     * @param colIndex column index
     * @return value of the matrix (plasma contribution only)
     */
    std::complex<double> getEntry(int rowIndex, int colIndex);

private:
    const int m_gridRes;
    const int m_NDOF;
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

};


#endif //HOTPLASMAKERNEL_MATRIX_H
