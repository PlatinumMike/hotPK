//
// Created by machiels on 29/08/2022.
//

#include "Mesh.h"

Mesh::Mesh(const int resolution, const double RWest, const double REast, const double RAnt) : m_res(resolution),
                                                                                              m_RWest(RWest),
                                                                                              m_REast(REast),
                                                                                              m_RAnt(RAnt),
                                                                                              m_elem(resolution - 1) {
    //todo: generate mesh, and add functions to convert indices local2global and so on.
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
