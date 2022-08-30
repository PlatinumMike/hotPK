#include "Engine.h"

int main() {

    //All inputs are in SI units! So also temperature (which is in energy units here) is in Joule, not eV.
    Engine engine;
    engine.run();

    return 0;
}
