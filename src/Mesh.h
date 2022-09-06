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

    /**
     * Retrieve list of element indices belonging to a given node.
     * @param nodeIndex index of the node
     * @param elemIndices vector of ints, will be filled
     */
    void getElemList(int nodeIndex, std::vector<int> &elemIndices);

    [[nodiscard]] double getNodePosition(int iGrid) const;
    [[nodiscard]] double getElemLeftPoint(int elemIndex) const;
    [[nodiscard]] double getElemRightPoint(int elemIndex) const;
    [[nodiscard]] double getElemMidPoint(int elemIndex) const;
    [[nodiscard]] double getElemWidth(int elemIndex) const;

    [[nodiscard]] bool isOnBoundary(int nodeIndex) const;

    int getElemCount();

    /**
     * Get "sector", so check if sample point R is left, inside, or right of given element.
     * @param R position R
     * @param elem element index
     * @return sector
     */
    int selectSector(double R, int elem);

    double getSValue(double R, int elem, double angle, int index);

    double getScaling(double R, int elem, int nodeIndex, double angle, int power);

    void saveNodes();

private:
    const int m_res; // number of nodes
    const int m_elem; //number of elements
    const double m_RWest;
    const double m_REast;
    const double m_RAnt;

    std::vector<double> m_nodePositions; //positions of all the mesh nodes
    std::vector<double> m_elementMids; //positions of all the element mid points
    std::vector<double> m_elementWidths; //positions of all the widths of the elements

    /**
     * get global node index from index local to the element
     * @param elemIndex index of the element
     * @param localNodeIndex left (0) or right (1)
     * @return global node
     */
    [[nodiscard]] int localNode2Global(int elemIndex, int localNodeIndex) const;

    /**
     * check if the selected element contains the left side of the basis function at node nodeIndex, or the right side
     * @param nodeIndex
     * @param elemIndex
     * @return true if the left side (so the rising part) is in the element.
     * @warning it assumes the left or right part is the element, so no check to see if it is fully outside the element.
     */
    bool isLeft(int nodeIndex, int elemIndex);

};


#endif //HOTPLASMAKERNEL_MESH_H
