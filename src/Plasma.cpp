//
// Created by machiels on 29/08/2022.
//

#include "Plasma.h"
#include "physicsConstants.h"
#include "AuxiliaryFunctions.h"

Plasma::Plasma(Mesh &mesh, int nTor, double omega, plasmaType pType) {
    m_mesh = &mesh;
    m_NSpecies = 0;
    m_nTor = nTor;
    m_omega = omega;
    m_pType = pType;
}

std::complex<double> Plasma::getCurrentMatrix(double R, int row, int col, int nodeIndex) {
    if (m_pType == hot) {
        return getCurrentMatrixHot(R, row, col, nodeIndex);
    }
    /*if the response is cold or warm, the response is local in space, so we use the 3x3 conductivity.
    So expressing plasma current J_p in terms of E,and expressing that in terms of \nabla Pot,A.
    we then take the derivative of the basis function explicitly. For hot plasma it is easier to
    compute the plasma response for just one, given basis function, so this uses the 3x4 conductivity,
    and no derivatives.*/
    if (col == 0) {
        std::complex<double> sigmaiR{0};
        std::complex<double> sigmaiphi{0};
        for (const auto &spec: specList) {
            sigmaiR += spec.getConductivity(R, row, 0);
            sigmaiphi += spec.getConductivity(R, row, 1);
        }
        return -physConstants::speedOfLight * (sigmaiR * m_mesh->tentDerivative(R, nodeIndex) +
                                               std::complex<double>{0, m_nTor / R} * sigmaiphi *
                                               m_mesh->tent(R, nodeIndex));
    } else {
        std::complex<double> sigmaij{0};
        for (const auto &spec: specList) {
            sigmaij += spec.getConductivity(R, row, col - 1);
        }
        return std::complex<double>{0, 1} * m_omega * sigmaij * m_mesh->tent(R, nodeIndex);
    }
}

void Plasma::addSpecies(double mass, double charge, double fraction, double omega, double peakTemp, plasmaType pType) {
    Species spec(mass, charge, fraction, m_nTor, omega, peakTemp, pType);
    specList.push_back(spec);
    m_NSpecies++;
}

std::complex<double> Plasma::getCurrentMatrixHot(double R, int row, int col, int nodeIndex) {
    std::complex<double> ans{0};
    constexpr int angularResolution = 101; //Make input later on.
    double angles[angularResolution];
    std::complex<double> integrand[angularResolution];
    std::vector<int> elemList{};
    m_mesh->getElemList(nodeIndex, elemList);

    for (int elem: elemList) {
        //get angle range
        int sector = m_mesh->selectSector(R, elem);
        AuxiliaryFunctions::getAngle(angles, angularResolution, sector);
        for (int angle = 0; angle < angularResolution; angle++) {
            integrand[angle] = 0; //zero first
            //get s1,s2
            double s1 = m_mesh->getSValue(R, elem, angles[angle], 0);
            double s2 = m_mesh->getSValue(R, elem, angles[angle], 1);
            for (int power = 1; power < 3; power++) {
                //sum over terms \propto s and \propto s^2.
                double scaling = m_mesh->getScaling(R, elem, nodeIndex, angles[angle], power);
                for (const auto &spec: specList) {
                    integrand[angle] += spec.getKernel(R, row, col, s1, s2, angles[angle], power) * scaling;
                }
            }
        }
        ans += AuxiliaryFunctions::integrateSimpson(angles, integrand, angularResolution);
    }

    return ans;
}