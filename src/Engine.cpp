//
// Created by machiels on 29/08/2022.
//

#include "Engine.h"
#include "Matrix.h"
#include <iostream>
#include "physicsConstants.h"
#include "omp.h"
#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"
//#include <filesystem> //not supported on my system

#ifndef GIT_COMMIT_HASH
#define GIT_COMMIT_HASH "?"
#endif

constexpr double pi = 3.141592653589793;

Engine::Engine(const std::string &fileName) {
    //initialize
    std::cout << "Git commit hash: " << GIT_COMMIT_HASH << std::endl;
#pragma omp parallel default(none)
    {
#pragma omp single nowait
        {
            std::printf("Starting program, using %d OpenMP threads\n", omp_get_num_threads());
        }
    }
    readInput(fileName);
}

void Engine::run() {
    const Parameters constInputs = inputs; //make a const copy, to avoid accidental modification of input variables.
    plasmaType pType;
    switch (constInputs.method) {
        case 0:
            pType = cold;
            std::printf("Using cold plasma model\n\n");
            break;
        case 1:
            pType = warm;
            std::printf("Using warm plasma model\n\n");
            break;
        default:
            pType = hot;
            std::printf("Using hot plasma model\n\n");
    }

    mesh_ptr = new Mesh(constInputs.gridResolution, constInputs.R_west, constInputs.R_east, constInputs.R_ant);
    plasma_ptr = new Plasma(*mesh_ptr, constInputs.nToroidal, constInputs.omega, constInputs.angularResolution, constInputs.R_east, constInputs.R_west, pType);

    double hydrogenDens = constInputs.electronDensity * constInputs.hydrogenConcentration;
    double hydrogenDensMin = constInputs.minElecDensity * constInputs.hydrogenConcentration;
    // add electrons
    plasma_ptr->addSpecies(physConstants::massElectron, -physConstants::elementaryCharge, constInputs.electronDensity,
                           constInputs.minElecDensity, constInputs.omega, constInputs.electronTemp, constInputs.minTemp,
                           constInputs.R0, constInputs.B0, pType);
    // add protons
    plasma_ptr->addSpecies(physConstants::massProton, physConstants::elementaryCharge, hydrogenDens, hydrogenDensMin,
                           constInputs.omega, constInputs.electronTemp, constInputs.minTemp, constInputs.R0,
                           constInputs.B0, pType);
    // add deuterons
    plasma_ptr->addSpecies(physConstants::massProton * 2, physConstants::elementaryCharge,
                           constInputs.electronDensity - hydrogenDens, constInputs.minElecDensity - hydrogenDensMin,
                           constInputs.omega, constInputs.electronTemp, constInputs.minTemp, constInputs.R0,
                           constInputs.B0, pType);

    matrix_ptr = new Matrix(constInputs.gridResolution, *mesh_ptr, *plasma_ptr, constInputs.nToroidal,
                            constInputs.omega, constInputs.J0R, constInputs.J0Phi, constInputs.J0Z);


    matrix_ptr->buildMatrix();

    if (constInputs.saveMatrix) {
        double startTime = omp_get_wtime();
        std::cout << "Writing data to disk now" << std::endl;
        matrix_ptr->saveMatrix("matRe.csv", true);
        matrix_ptr->saveMatrix("matIm.csv", false);
        double endTime = omp_get_wtime();
        std::cout << "Finished writing, time used: " << endTime - startTime << std::endl;
    }

    if (constInputs.solveMatrix) {
        matrix_ptr->solve();
        matrix_ptr->saveSolution("sol.csv");
        mesh_ptr->saveNodes();
    }
}

int Engine::readInput(const std::string &inputFileName) {
    std::cout << "Reading input from file " << inputFileName << std::endl;
    if (true) { //std::filesystem::exists(inputFileName)) {
        boost::property_tree::ptree root;
        boost::property_tree::read_json(inputFileName, root);

        //todo: split up struct into multiple pieces, not all these variables need to be known in every class.
        //create new instance of Parameters, then fill it up using the json file.
        inputs.gridResolution = root.get<int>("resolution");
        inputs.R_west = root.get<double>("R_west");
        inputs.R_east = root.get<double>("R_east");
        inputs.R0 = root.get<double>("R0");
        inputs.R_ant = root.get<double>("R_ant");

        inputs.nToroidal = root.get<int>("nToroidal");
        inputs.B0 = root.get<double>("B0");
        inputs.omega = root.get<double>("omega");
        inputs.J0R = root.get<double>("J0R");
        inputs.J0Phi = root.get<double>("J0Phi");
        inputs.J0Z = root.get<double>("J0Z");

        inputs.saveMatrix = root.get<bool>("saveMatrix");
        inputs.solveMatrix = root.get<bool>("solveMatrix");

        inputs.electronDensity = root.get<double>("electronDensity");
        inputs.minElecDensity = root.get<double>("minElecDensity");
        inputs.hydrogenConcentration = root.get<double>("hydrogenConcentration");
        inputs.electronTemp = root.get<double>("electronTemp");
        inputs.minTemp = root.get<double>("minTemp");
        inputs.angularResolution = root.get<int>("angularResolution");
        inputs.method = root.get<int>("method");

        //compute derived quantities
        inputs.R_center = 0.5 * (inputs.R_west + inputs.R_east);
        inputs.numDegFreedom = inputs.gridResolution * 4;
        inputs.antennaPeriod = 2.0 * pi / inputs.omega;
        inputs.omegaOverC = inputs.omega / physConstants::speedOfLight;

        std::printf("Grid resolution: %d\n", inputs.gridResolution);
        std::printf("Number of degrees of freedom: %d\n", inputs.numDegFreedom);

        return 0;
    } else {
        std::cout << "Cannot find input file " << inputFileName << std::endl;
        exit(EXIT_FAILURE);
        return 1;
    }
}
