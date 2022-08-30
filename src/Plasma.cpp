//
// Created by machiels on 29/08/2022.
//

#include "Plasma.h"
#include "physicsConstants.h"

Plasma::Plasma(Mesh &mesh, double nTor, double omega) {
    m_mesh = &mesh;
    m_NSpecies = 0;
    m_nTor = nTor;
    m_omega = omega;
}

std::complex<double> Plasma::getCurrentMatrix(double R, int row, int col, int nodeIndex) {
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

void Plasma::addSpecies(double mass, double charge, double fraction, double omega) {
    Species spec(mass, charge, fraction, omega);
    specList.push_back(spec);
    m_NSpecies++;
}
