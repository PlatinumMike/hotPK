cmake_minimum_required(VERSION 3.21)
project(hotPlasmaKernel)

set(CMAKE_CXX_STANDARD 17)

set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}")

add_executable(hotPlasmaKernel src/main.cpp src/Mesh.cpp src/Mesh.h src/Matrix.cpp src/Matrix.h src/Plasma.cpp src/Plasma.h src/Engine.cpp src/Engine.h)

message("Build type       : " ${CMAKE_BUILD_TYPE}) #Prints Debug or Release, so it's clear if the code has actually been compiled with/without optimizations
message("Eigen dir        : " ${EIGEN3_INCLUDE_DIR}) #there actually is a more modern way using find_package but that did not work for me