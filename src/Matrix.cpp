//
// Created by machiels on 29/08/2022.
//

#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"
#include "Writers.h"
#include <boost/math/quadrature/gauss.hpp>
#include "omp.h"

Matrix::Matrix(int gridRes, Mesh &mesh, Plasma &plasma, int nTor, double omega) : m_gridRes(gridRes), m_NDOF(4 * gridRes),
                                                                    m_nTor(nTor), m_omega(omega){
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
    double startTime = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp for collapse(2)
        for (int i = 0; i < m_NDOF; ++i) {
            for (int j = 0; j < m_NDOF; ++j) {
                (*globalMatrix)(i, j) = getEntry(i, j);
            }
        }
#pragma omp for
        for (int i = 0; i < m_NDOF; ++i) {
            (*rhs)(i) = getRhs(i);
        }
    }
    double endTime = omp_get_wtime();
    std::cout << "Matrix build completed, time used: " << endTime - startTime << std::endl;
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
    int iNode = global2Node(rowIndex);
    int jNode = global2Node(colIndex);
    int iComp = global2Comp(rowIndex);
    int jComp = global2Comp(colIndex);

    // boundary check
    if (iNode == 0 || iNode == m_gridRes - 1) {
        // on the boundary, so for rows associated with Pot, Aphi or AZ, return 0. This is because those will be overruled by essential BCs.
        if (iComp != 1) {
            return {0};
        }
    }
    // in any case, there is never any dependency on Pot, Aphi, AZ that lie on the boundary, as they are 0. So also zero out columns related to those.
    if (jNode == 0 || jNode == m_gridRes - 1) {
        if (jComp != 1) {
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