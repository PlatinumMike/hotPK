//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MESH_H
#define HOTPLASMAKERNEL_MESH_H

#include <vector>


class Mesh {
public:
    Mesh(int resolution, double RWest, double REast, double RAnt);
    /**
     * Compute basis function
     * @param R major radius
     * @param nodeIndex index of the grid node
     * @return return value of basis function at position R
     */
    double tent(double R, int nodeIndex);
    /**
     * d/dR of basis function
     */
    double tentDerivative(double R, int nodeIndex);

    double getRAnt();

private:
    const int m_res; // number of nodes
    const int m_elem; //number of elements
    const double m_RWest;
    const double m_REast;
    const double m_RAnt;

    std::vector<double> m_nodePositions; //positions of all the mesh nodes
    std::vector<double> m_elementMids; //positions of all the element mid points
    std::vector<double> m_elementWidths; //positions of all the widths of the elements

    [[nodiscard]] double getNodePosition(int iGrid) const;
    [[nodiscard]] double getElemMidPoint(int elemIndex) const;
    [[nodiscard]] double getElemWidth(int elemIndex) const;

};


#endif //HOTPLASMAKERNEL_MESH_H
