//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MESH_H
#define HOTPLASMAKERNEL_MESH_H

#include <vector>


class Mesh {
public:
    /**
     * Initialize mesh object
     * @param resolution mesh resolution
     * @param RWest west wall position (major radius)
     * @param REast east wall position
     * @param RAnt antenna position
     */
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

    /**
     * get antenna position (just a single current sheet)
     * @return major radius of antenna
     */
    double getRAnt();

    /**
     * Retrieve list of element indices belonging to a given node.
     * @param nodeIndex index of the node
     * @param elemIndices vector of ints, will be filled
     */
    void getElemList(int nodeIndex, std::vector<int> &elemIndices);

    /**
     * get position of a mesh node
     * @param iGrid node index
     * @return node position
     */
    [[nodiscard]] double getNodePosition(int iGrid) const;
    /**
     * get left wall position of element
     * @param elemIndex index of element
     * @return element left position
     */
    [[nodiscard]] double getElemLeftPoint(int elemIndex) const;
    /**
     * get right wall position of element
     * @param elemIndex index of element
     * @return element right position
     */
    [[nodiscard]] double getElemRightPoint(int elemIndex) const;
    /**
     * get geometric middle position of element
     * @param elemIndex index of element
     * @return element mid position
     */
    [[nodiscard]] double getElemMidPoint(int elemIndex) const;
    /**
     * get width of element
     * @param elemIndex index of element
     * @return element width
     */
    [[nodiscard]] double getElemWidth(int elemIndex) const;

    /**
     * check if node is on the metal wall (so on the domain boundary)
     * @param nodeIndex index of the mesh node
     * @return true if it is on the boundary, false otherwise
     */
    [[nodiscard]] bool isOnBoundary(int nodeIndex) const;

    /**
     * @return total number of elements
     */
    int getElemCount();

    /**
     * Get "sector", so check if sample point R is left, inside, or right of given element.
     * @param R position R
     * @param elem element index
     * @return sector
     */
    int selectSector(double R, int elem);

    /**
     * Get integral bound in radial direction, called "s"
     * @param R major radius
     * @param elem element index
     * @param angle polar angle
     * @param index select \f$ s_1, s_2, s_3, s_4 \f$
     * @return \f$ s \f$
     */
    double getSValue(double R, int elem, double angle, int index);

    /**
     * computes extra prefactor related to the basis functions
     * @param R major radius
     * @param elem element index
     * @param nodeIndex index of the mesh node
     * @param angle polar angle
     * @param power toggle between the two contributions, which have a different power of "s" in the integrand
     * @return prefactor
     */
    double getScaling(double R, int elem, int nodeIndex, double angle, int power);

    /**
     * write mesh data to disk
     */
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
