//
// Created by machiels on 29/08/2022.
//

#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"

constexpr double pi = 3.141592653589793;

Matrix::Matrix(int gridRes, Mesh &mesh, int nTor, double antFreq) : m_gridRes(gridRes), m_NDOF(4 * gridRes),
                                                                    m_nTor(nTor), m_omega(2 * pi * antFreq) {
    globalMatrix = new Eigen::MatrixXcd(m_NDOF,m_NDOF);
    rhs = new Eigen::VectorXcd(m_NDOF);
    solution = new Eigen::VectorXcd(m_NDOF);
    m_mesh = &mesh;
}

std::complex<double> Matrix::getEntry(int rowIndex, int colIndex) {
    return {1};
}

int Matrix::global2Node(int i) {
    return i%m_gridRes;
}

int Matrix::global2Comp(int i) {
    //components are ordered as such: Pot, Ar, Aphi, Az
    return i/m_gridRes;
}

void Matrix::buildMatrix() {
    //todo: use OMP to multi-thread?
    for (int i = 0; i < m_NDOF; ++i) {
        for (int j = 0; j < m_NDOF; ++j) {
            (*globalMatrix)(i, j) = getEntry(i, j);
        }
    }
    for (int i = 0; i < m_NDOF; ++i) {
        (*rhs)(i) = getRhs(i);
    }
}

void Matrix::solve() {
    //LU factorization, in place.
    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXcd> > lu(*globalMatrix);
    *solution = lu.solve(*rhs);
    std::cout<<*solution<<std::endl;
}

std::complex<double> Matrix::getRhs(int rowIndex) {
    int node = global2Node(rowIndex);
    int comp = global2Comp(rowIndex);
    //antenna current. Todo: make this an input
    constexpr double JantR = 0.0;
    constexpr double JantPhi = 0.0;
    constexpr double JantZ = 1.0;
    double Rant = m_mesh->getRAnt();
    if (comp == 0) {
        //Potential
        return -std::complex<double>{0, physConstants::mu_0 * physConstants::speedOfLight / m_omega} * (
                m_mesh->tentDerivative(Rant, node) * JantR -
                std::complex<double>{0, m_nTor / Rant} * m_mesh->tent(Rant, node) * JantPhi);
    } else if (comp == 1) {
        //AR
        return -physConstants::mu_0 * JantR * m_mesh->tent(Rant, node);
    } else if (comp == 2) {
        //Aphi
        return -physConstants::mu_0 * JantPhi * m_mesh->tent(Rant, node);
    } else {
        //AZ
        return -physConstants::mu_0 * JantZ * m_mesh->tent(Rant, node);
    }
}
Matrix::~Matrix(){
    delete globalMatrix;
    delete rhs;
    delete solution;
    delete m_mesh;
}
