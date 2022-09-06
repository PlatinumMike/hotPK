//
// Created by machiels on 29/08/2022.
//

#include "Engine.h"
#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"
#include "omp.h"

constexpr double pi = 3.141592653589793;

Engine::Engine() {
    constexpr double antennaFreq = 40e6;
    double omega = 2 * pi * antennaFreq;
    int nTor = 27;
    double RWest = 2.3;
    double REast = 4.0;
    double RAnt = 3.9;
    int resolution = 300;
    plasmaType pType = hot;
    double peakTemp = 4.0e3*physConstants::elementaryCharge;

    std::cout<< "Git commit hash: "<<GIT_COMMIT_HASH<<std::endl;
#pragma omp parallel default(none)
    {
#pragma omp single nowait
        {
            std::printf("Starting program, using %d OpenMP threads\n",omp_get_num_threads());
        }
    }

    mesh_ptr = new Mesh(resolution, RWest, REast, RAnt);
    plasma_ptr = new Plasma(*mesh_ptr, nTor, omega, pType);
    plasma_ptr->addSpecies(physConstants::massElectron, -physConstants::elementaryCharge, 1.0, omega, peakTemp, pType);
    plasma_ptr->addSpecies(physConstants::massProton, physConstants::elementaryCharge, 1.0e-2, omega, peakTemp, pType);
    plasma_ptr->addSpecies(physConstants::massProton * 3, physConstants::elementaryCharge * 2, 0.495, omega, peakTemp, pType);
    matrix_ptr = new Matrix(resolution, *mesh_ptr, *plasma_ptr, nTor, omega);
}

void Engine::run() {
    matrix_ptr->buildMatrix();

    double startTime = omp_get_wtime();
    std::cout<<"Writing data to disk now"<<std::endl;
    matrix_ptr->saveMatrix("matRe.csv", true);
    matrix_ptr->saveMatrix("matIm.csv", false);
    double endTime = omp_get_wtime();
    std::cout<<"Finished writing, time used: "<< endTime - startTime << std::endl;

    matrix_ptr->solve();
    matrix_ptr->saveSolution("sol.csv");
    mesh_ptr->saveNodes();
}
