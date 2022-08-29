cmake_minimum_required(VERSION 3.21)
project(hotPlasmaKernel)

set(CMAKE_CXX_STANDARD 17)

add_executable(hotPlasmaKernel src/main.cpp src/Mesh.cpp src/Mesh.h src/Matrix.cpp src/Matrix.h src/Plasma.cpp src/Plasma.h src/Engine.cpp src/Engine.h)
