//
// Created by machiels on 29/08/2022.
//

#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"
#include "Writers.h"
#include <boost/math/quadrature/gauss.hpp>

Matrix::Matrix(int gridRes, Mesh &mesh, Plasma &plasma, int nTor, double omega, plasmaType pType) : m_gridRes(gridRes), m_NDOF(4 * gridRes),
                                                                    m_nTor(nTor), m_omega(omega), m_pType(pType) {
    globalMatrix = new Eigen::MatrixXcd(m_NDOF,m_NDOF);
    rhs = new Eigen::VectorXcd(m_NDOF);
    solution = new Eigen::VectorXcd(m_NDOF);
    m_mesh = &mesh;
    m_plasma = &plasma;
}

std::complex<double> Matrix::getEntry(int rowIndex, int colIndex) {
    return getEntryVacuum(rowIndex, colIndex) + getEntryPlasma(rowIndex, colIndex);
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
}

std::complex<double> Matrix::getEntryPlasma(int rowIndex, int colIndex){
    if (m_pType == vacuum){
        return {0};
    }
    //else:

    int iNode = global2Node(rowIndex);
    int jNode = global2Node(colIndex);
    int iComp = global2Comp(rowIndex);
    int jComp = global2Comp(colIndex);

    // boundary check
    if (iNode == 0 || jNode == 0 || iNode == m_gridRes-1 || jNode == m_gridRes -1){
        // if on boundary, zero contribution to Pot, Aphi, Az. But not that of Ar.
        if (iComp != 1){
            return {0};
        }
    }

    std::complex<double> returnVal{0};
    auto func = [iNode, jNode, iComp, jComp, this](double R) {
        if (iComp == 0) {
            return physConstants::mu_0 * R * std::complex<double>{0, 1} * physConstants::speedOfLight / m_omega *
                   (m_mesh->tentDerivative(R, iNode) * m_plasma->getCurrentMatrix(R, 0, jComp, jNode) -
                    std::complex<double>{0, m_nTor / R} * m_mesh->tent(R, iNode) *
                    m_plasma->getCurrentMatrix(R, 1, jComp, jNode));
        } else {
            return physConstants::mu_0 * R * m_mesh->tent(R, iNode) *
                   m_plasma->getCurrentMatrix(R, iComp - 1, jComp, jNode);
        }
    };

    /*if (m_pType == cold){
        return {0};
    }else if (m_pType == warm){
        return {0}; //todo
    }else{
        //hot plasma, todo
        return {0};
    }*/

    std::vector<int> elemIndicesI{};
    m_mesh->getElemList(iNode,elemIndicesI);

    for(int const &elemI : elemIndicesI){
        double Rleft = m_mesh->getElemLeftPoint(elemI);
        double Rright = m_mesh->getElemRightPoint(elemI);
        //5-point Gauss-Legendre quadrature
        returnVal += boost::math::quadrature::gauss<double,5>::integrate(func,Rleft,Rright);
    }
    return returnVal;
}

std::complex<double> Matrix::getEntryVacuum(int rowIndex, int colIndex) {
    return {0}; //todo: placeholder. Should return the vacuum contribution.
}

void Matrix::saveMatrix(std::string fileName, bool isReal) {
    Writers::writeCsv2(fileName,*globalMatrix, isReal);
}
