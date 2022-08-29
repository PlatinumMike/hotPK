//
// Created by machiels on 29/08/2022.
//

#include "Engine.h"
#include "Matrix.h"
#include <iostream>

Engine::Engine(){
    matrix_ptr = new Matrix(30);
}

void Engine::run() {
    //todo: for now some random indices
    std::cout<< matrix_ptr->getEntry(10,4) << std::endl;
}
