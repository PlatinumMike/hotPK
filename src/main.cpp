#include "Engine.h"
#include <iostream>

int main(int argc, char *argv[]) {
    //All inputs are in SI units! So also temperature (which is in energy units here) is in Joule, not eV.
    if (argc == 1) {
        std::cout << "Insufficient command line arguments, json input file name missing!" << std::endl;
        return 1;
    } else if (argc > 2) {
        std::cout << "Too many command line arguments!" << std::endl;
        return 2;
    } else {
        std::string inputFileName(argv[1]);
        Engine engine(inputFileName); //initialize engine
        engine.run(); //start the solving process
        return 0;
    }
}
