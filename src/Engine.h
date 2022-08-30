//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_ENGINE_H
#define HOTPLASMAKERNEL_ENGINE_H

#include "Matrix.h"
#include "Mesh.h"
#include "Plasma.h"

class Engine {
public:
    Engine();
    void run();

private:
    Mesh *mesh_ptr;
    Plasma *plasma_ptr;
    Matrix *matrix_ptr;

};


#endif //HOTPLASMAKERNEL_ENGINE_H
