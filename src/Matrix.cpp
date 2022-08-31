//
// Created by machiels on 29/08/2022.
//

#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"
#include "Writers.h"
#include <boost/math/quadrature/gauss.hpp>
#include "omp.h"
#include "ElementIntegrals.h"

Matrix::Matrix(int gridRes, Mesh &mesh, Plasma &plasma, int nTor, double omega) : m_gridRes(gridRes), m_NDOF(4 * gridRes),
        m_nTor(nTor), m_omega(omega), m_omegaOverC2(omega*omega/(physConstants::speedOfLight*physConstants::speedOfLight)){
    globalMatrix = new Eigen::MatrixXcd(m_NDOF,m_NDOF);
    rhs = new Eigen::VectorXcd(m_NDOF);
    solution = new Eigen::VectorXcd(m_NDOF);
    m_mesh = &mesh;
    m_plasma = &plasma;
    sparseMatrix = new Eigen::SparseMatrix<std::complex<double>>(m_NDOF,m_NDOF);

    buildRHS();
    buildSparse();
}

std::complex<double> Matrix::getEntry(int rowIndex, int colIndex) {
    return getEntryVacuum(rowIndex, colIndex) + getEntryPlasma(rowIndex, colIndex);
}

int Matrix::global2Node(int i) const {
    return i%m_gridRes;
}

int Matrix::global2Comp(int i) const {
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
    }
    double endTime = omp_get_wtime();
    std::cout << "Matrix build completed, time used: " << endTime - startTime << std::endl;
}

void Matrix::solve() {
    double startTime = omp_get_wtime();
    //LU factorization, in place.
    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXcd> > lu(*globalMatrix);
    *solution = lu.solve(*rhs);
    double endTime = omp_get_wtime();
    std::cout << "Solve completed, time used: " << endTime - startTime << std::endl;
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
    if (m_mesh->isOnBoundary(iNode)) {
        // on the boundary, so for rows associated with Pot, Aphi or AZ, return 0. This is because those will be overruled by essential BCs.
        if (iComp != 1) {
            return {0};
        }
    }
    // in any case, there is never any dependency on Pot, Aphi, AZ that lie on the boundary, as they are 0. So also zero out columns related to those.
    if (m_mesh->isOnBoundary(jNode)) {
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
    return sparseMatrix->coeff(rowIndex,colIndex);
}

void Matrix::saveMatrix(const std::string& fileName, bool isReal) {
    Writers::writeCsv2(fileName,*globalMatrix, isReal);
}

void Matrix::buildRHS() {
    double startTime = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < m_NDOF; ++i) {
            (*rhs)(i) = getRhs(i);
        }
    }
    double endTime = omp_get_wtime();
    std::cout << "RHS build completed, time used: " << endTime - startTime << std::endl;
}

void Matrix::buildSparse() {
    double startTime = omp_get_wtime();
    std::vector<T> tripletList;
    int numElem = m_mesh->getElemCount();
    tripletList.reserve(numElem * 8 * 8); //estimate of nnz: for every element, a 8x8 block is inserted. (overlap is possible, but the setFromTriplets sums over duplicate indices).
    Eigen::Matrix<std::complex<double>,8,8>  localMatrix; //8x8 matrix block, consisting of the total contribution of one element to the global matrix. it is a 2x2 block of 4x4 subblocks. Each 4x4 subblock belongs to one node.
    Eigen::Matrix2cd miniMatrix; //the mini matrix is 2 by 2 for linear elements in 1D, 3x3 for linear triangles, nxn for a linear element with n corners.
    //element-wise construction, which is easier than node-wise in a way. Especially in 2D or 3D where you do not need to worry about the number of neighbours and so on.

    std::complex<double> imagUnit{0,1};
    for (int elem = 0; elem < numElem; ++elem) {
        localMatrix.setZero(); //zero each iteration
        ElementIntegrals integrals(m_mesh->getElemMidPoint(elem), m_mesh->getElemWidth(elem)); //get integrals

        //G component of wave equation
        miniMatrix = -m_nTor * m_nTor * integrals.moment1m1 - integrals.moment4p1; //GPot
        commitMiniMatrix(localMatrix, miniMatrix, 0, 0);
        miniMatrix = 0.5 * imagUnit * m_omega / physConstants::speedOfLight *
                     (integrals.moment2p1 - integrals.moment3p1 - integrals.moment1p0); //GAR
        commitMiniMatrix(localMatrix, miniMatrix, 0, 1);
        miniMatrix = m_nTor * m_omega / physConstants::speedOfLight * integrals.moment1p0; //GAphi
        commitMiniMatrix(localMatrix, miniMatrix, 0, 2);
        //GAZ, skip cause 0.

        //FR component of wave equation
        miniMatrix = 0.5 * imagUnit * m_omega / physConstants::speedOfLight *
                     (integrals.moment3p1 - integrals.moment2p1 - integrals.moment1p0); //FRPot
        commitMiniMatrix(localMatrix, miniMatrix, 1, 0);
        miniMatrix = -(m_nTor * m_nTor + 1) * integrals.moment1m1 + m_omegaOverC2 * integrals.moment1p1 -
                     integrals.moment2p0 - integrals.moment3p0 - integrals.moment4p1; //FRAR
        commitMiniMatrix(localMatrix, miniMatrix, 1, 1);
        miniMatrix =
                -m_nTor * (2 * integrals.moment1m1 + integrals.moment2p0 + integrals.moment3p0)* imagUnit; //FRAphi
        commitMiniMatrix(localMatrix, miniMatrix, 1, 2);
        //FRAZ, skip cause 0.

        //Fphi component of wave equation
        miniMatrix = -m_nTor * m_omega / physConstants::speedOfLight * integrals.moment1p0; //FphiPot
        commitMiniMatrix(localMatrix, miniMatrix, 2, 0);
        miniMatrix = m_nTor * (2 * integrals.moment1m1 + integrals.moment2p0 + integrals.moment3p0) * imagUnit; //FphiAR
        commitMiniMatrix(localMatrix, miniMatrix, 2, 1);
        miniMatrix = -(m_nTor * m_nTor + 1) * integrals.moment1m1 + m_omegaOverC2 * integrals.moment1p1 -
                     integrals.moment2p0 - integrals.moment3p0 - integrals.moment4p1; //FphiAphi
        commitMiniMatrix(localMatrix, miniMatrix, 2, 2);
        //FphiAZ, skip cause 0.


        //FZ component of wave equation
        //FZPot, skip cause 0.
        //FZAR, skip cause 0.
        //FZAphi, skip cause 0.
        miniMatrix = -m_nTor * m_nTor * integrals.moment1m1 + m_omegaOverC2 * integrals.moment1p1 -
                     integrals.moment4p1; //FZAZ
        commitMiniMatrix(localMatrix, miniMatrix, 3, 3);

        //now a few exceptions regarding the BCs
        if(elem==0){
            //west boundary element, skipping over the AR component, because that is a natural boundary condition, so this is already handled by the equation itself, so no modification needed, can be inserted as normal.
            //zero out rows belonging to Pot, Aphi, AZ, then adding essential Dirichlet BCs manually:
            for (int i = 0; i < 4; ++i) {
                if (i != 1) {
                    localMatrix.row(i).setZero(); //here you zero a row, but you can also zero the column, because we have homogeneous Dirichlet BCs, so AR will not depend on the other components on the boundary.
                    localMatrix.col(i).setZero();
                    localMatrix(i,i) = 1.0;
                    (*rhs)(local2Global(elem, 0, i)) = 0;
                }
            }

        }else if(elem==numElem-1){
            //east boundary element
            for (int i = 4; i < 8; ++i) {
                if (i != 5) {
                    localMatrix.row(i).setZero();
                    localMatrix.col(i).setZero();
                    localMatrix(i,i) = 1.0;
                    (*rhs)(local2Global(elem, 1, i-4)) = 0;
                }
            }
        }

        //ready to commit local matrix to the global one.
        commitLocalMatrix(tripletList, localMatrix, elem);

    }
    sparseMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
    double endTime = omp_get_wtime();
    std::cout << "Sparse matrix build completed, time used: " << endTime - startTime << std::endl;
}

void Matrix::commitMiniMatrix(Eigen::Matrix<std::complex<double>, 8, 8> &localMatrix,
                              const Eigen::Matrix2cd &miniMatrix, int comp1, int comp2) {
    //loop over both basis functions of the element.
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            localMatrix(comp1 + 4 * i, comp2 + 4 * j) += miniMatrix(i, j);
        }
    }
}

void Matrix::commitLocalMatrix(std::vector<T> &tripletList, const Eigen::Matrix<std::complex<double>, 8, 8> &localMatrix,
                          int elemIndex) {
    int row, col;
    //loop over all 4 subblocks, related to the basis functions of the element.
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            //loop over all entries within one subblock.
            for (int comp1 = 0; comp1 < 4; ++comp1) {
                for (int comp2 = 0; comp2 < 4; ++comp2) {
                    row = local2Global(elemIndex, i, comp1);
                    col = local2Global(elemIndex, j, comp2);
                    tripletList.push_back(T(row, col, localMatrix(comp1 + 4 * i, comp2 + 4 * j)));
                }
            }
        }
    }

}

int Matrix::local2Global(int elemIndex, int localNodeIndex, int component) {
    return component * m_gridRes + elemIndex + localNodeIndex;
}

void Matrix::saveSolution(const std::string& fileName) {
    Writers::writeCsv(fileName, *solution);
}
