//
// Created by machiels on 29/08/2022.
//

#include "Mesh.h"
#include <cmath>
#include "Writers.h"

constexpr double pi = 3.141592653589793;

Mesh::Mesh(const int resolution, const double RWest, const double REast, const double RAnt) : m_res(resolution),
                                                                                              m_RWest(RWest),
                                                                                              m_REast(REast),
                                                                                              m_RAnt(RAnt),
                                                                                              m_elem(resolution - 1) {
    m_nodePositions.resize(m_res);
    m_elementMids.resize(m_elem);
    m_elementWidths.resize(m_elem);

    //uniform grid spacing for now
    double dR = (m_REast - m_RWest) / m_elem;
    for (int i = 0; i < m_res; ++i) {
        m_nodePositions[i] = m_RWest + dR * i;
    }
    for (int i = 0; i < m_elem; ++i) {
        m_elementMids[i] = 0.5 * (m_nodePositions[i] + m_nodePositions[i + 1]);
        m_elementWidths[i] = m_nodePositions[i + 1] - m_nodePositions[i];
    }
}

double Mesh::getNodePosition(int iGrid) const {
    return m_nodePositions[iGrid];
}

double Mesh::getElemMidPoint(int elemIndex) const {
    return m_elementMids[elemIndex];
}

double Mesh::getElemWidth(int elemIndex) const {
    return m_elementWidths[elemIndex];
}

double Mesh::tent(double R, int nodeIndex) {
    if (R >= m_nodePositions[nodeIndex - 1] && R < m_nodePositions[nodeIndex]) {
        return (R - m_nodePositions[nodeIndex - 1]) / (m_nodePositions[nodeIndex] - m_nodePositions[nodeIndex - 1]);
    } else if (R >= m_nodePositions[nodeIndex] && R <= m_nodePositions[nodeIndex + 1]) {
        return (m_nodePositions[nodeIndex + 1] - R) / (m_nodePositions[nodeIndex + 1] - m_nodePositions[nodeIndex]);
    } else {
        return 0;
    }
}

double Mesh::tentDerivative(double R, int nodeIndex) {
    if (R >= m_nodePositions[nodeIndex - 1] && R < m_nodePositions[nodeIndex]) {
        return 1.0 / (m_nodePositions[nodeIndex] - m_nodePositions[nodeIndex - 1]);
    } else if (R >= m_nodePositions[nodeIndex] && R <= m_nodePositions[nodeIndex + 1]) {
        return -1.0 / (m_nodePositions[nodeIndex + 1] - m_nodePositions[nodeIndex]);
    } else {
        return 0;
    }
}

double Mesh::getRAnt() {
    return m_RAnt;
}

void Mesh::getElemList(int nodeIndex, std::vector<int> &elemIndices) {
    if(nodeIndex==0){
        elemIndices.push_back(0);
    }else if (nodeIndex==m_res-1){
        elemIndices.push_back(m_elem-1);
    }else{
        elemIndices.push_back(nodeIndex-1);
        elemIndices.push_back(nodeIndex);
    }
}

double Mesh::getElemLeftPoint(int elemIndex) const {
    return m_nodePositions[localNode2Global(elemIndex,0)];
}

double Mesh::getElemRightPoint(int elemIndex) const {
    return m_nodePositions[localNode2Global(elemIndex,1)];
}

int Mesh::localNode2Global(int elemIndex, int localNodeIndex) const {
    return elemIndex + localNodeIndex; //assuming the mesh is ordered.
}

bool Mesh::isOnBoundary(int nodeIndex) const {
    return (nodeIndex == 0 || nodeIndex == m_res - 1);
}

int Mesh::getElemCount() {
    return m_elem;
}

int Mesh::selectSector(double R, int elem) {
    double Rleft = getElemLeftPoint(elem);
    double Rright = getElemRightPoint(elem);
    if (R <= Rleft) {
        return 1; //left
    } else if (R < Rright) {
        return 2; //middle
    } else {
        return 3; //right
    }
}

double Mesh::getSValue(double R, int elem, double angle, int index) {
    double sValue;
    double Rleft = getElemLeftPoint(elem);
    double Rright = getElemRightPoint(elem);
    double sin = std::sin(angle);
    constexpr double lowerLimit = 1.0e-12;
    if (std::abs(sin) < lowerLimit) {
        sin = lowerLimit; //in order to avoid div by 0. if \sin(\alpha) \to 0, s \to \infty, but if it is 10^30, or 10^300 is pretty much irrelevant cause the H functions will have vanished long before that. So justified to have a lower limit.
    }
    if (index == 0) {
        //s1 (or s3)
        if (R <= Rleft) {
            sValue = (R - Rleft) / sin;
        } else if (R < Rright) {
            sValue = 0;
        } else {
            sValue = (R - Rright) / sin;
        }
    } else {
        //s2 (or s4)
        if (R <= Rleft) {
            sValue = (R - Rright) / sin;
        } else if (R < Rright) {
            if (angle > 0) {
                sValue = (R - Rleft) / sin;
            } else {
                sValue = (R - Rright) / sin;
            }
        } else {
            sValue = (R - Rleft) / sin;
        }
    }
    return std::abs(sValue); //need to take abs() in the edge case where sin is close to 0. This is because then the sin is set to lower limit so the information of its sign is lost. In any case, we know that s>=0, so free to take abs().
}

bool Mesh::isLeft(int nodeIndex, int elemIndex) {
    if (nodeIndex == elemIndex + 1){
        return true;
    }else if (nodeIndex == elemIndex){
        return false;
    }else{
        return false; //todo: placeholder, this should actually throw an error because the check failed, element contains no part of basis at nodeIndex.
    }
}

double Mesh::getScaling(double R, int elem, int nodeIndex, double angle, int power) {
    double Rleft = getElemLeftPoint(elem);
    double Rright = getElemRightPoint(elem);
    double Rj = getNodePosition(nodeIndex);
    if (isLeft(nodeIndex, elem)) {
        if (power == 1) {
            return (R - Rleft) / (Rj - Rleft);
        } else {
            return -std::sin(angle) / (Rj - Rleft);
        }
    } else {
        if (power == 1) {
            return (Rright - R) / (Rright - Rj);
        } else {
            return std::sin(angle) / (Rright - Rj);
        }
    }
}

void Mesh::saveNodes() {
    Writers::writeCsv("nodes.csv", m_nodePositions.data(), m_res);
}
