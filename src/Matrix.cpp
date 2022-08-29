//
// Created by machiels on 29/08/2022.
//

#include "Matrix.h"

Matrix::Matrix(int gridRes) : m_gridRes(gridRes), m_NDOF(4*gridRes) {

}

std::complex<double> Matrix::getEntry(int rowIndex, int colIndex) {
    return {1};
}

int Matrix::global2Node(int i) {
    return i%m_gridRes;
}

int Matrix::global2Comp(int i) {
    return i/m_gridRes;
}
