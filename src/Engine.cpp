//
// Created by machiels on 29/08/2022.
//

#include "Engine.h"
#include "Matrix.h"
#include <iostream>

Engine::Engine(){
    mesh_ptr = new Mesh(30,2.3,4.0,3.9);
    matrix_ptr = new Matrix(30,*mesh_ptr,27,40e6);
}

void Engine::run() {
    //todo: for now some random indices
    std::cout<< matrix_ptr->getEntry(10,4) << std::endl;
    matrix_ptr->buildMatrix();
    //matrix_ptr->solve();
}
