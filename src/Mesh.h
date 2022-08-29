//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_MESH_H
#define HOTPLASMAKERNEL_MESH_H


class Mesh {
public:
    Mesh(int resolution, double RWest, double REast, double RAnt);

private:
    int m_res;
    double m_RWest;
    double m_REast;
    double m_RAnt;

};


#endif //HOTPLASMAKERNEL_MESH_H
