//
// Created by machiels on 29/08/2022.
//

#ifndef HOTPLASMAKERNEL_ENGINE_H
#define HOTPLASMAKERNEL_ENGINE_H

#include "Matrix.h"
#include "Mesh.h"
#include "Plasma.h"
#include "Parameters.h"
#include <string>

class Engine {
public:
    explicit Engine(const std::string &fileName);
    void run();

private:
    Mesh *mesh_ptr;
    Plasma *plasma_ptr;
    Matrix *matrix_ptr;
    Parameters inputs;
    int readInput(const std::string &inputFileName);

};


#endif //HOTPLASMAKERNEL_ENGINE_H
