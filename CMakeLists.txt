cmake_minimum_required(VERSION 3.21)
project(hotPlasmaKernel)

set(CMAKE_CXX_STANDARD 17)

set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

set( BOOST_INCLUDE_DIR "$ENV{BOOST_INCLUDE_DIR}" )
IF( NOT BOOST_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable BOOST_INCLUDE_DIR to the include directory of your Boost installation.")
ENDIF()

INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" "${BOOST_INCLUDE_DIR}")

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

add_executable("${PROJECT_NAME}.x" src/main.cpp src/Mesh.cpp src/Mesh.h src/Matrix.cpp src/Matrix.h src/Plasma.cpp src/Plasma.h src/Engine.cpp src/Engine.h src/Writers.cpp src/Writers.h src/Species.cpp src/Species.h)

target_compile_definitions("${PROJECT_NAME}.x" PRIVATE "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries("${PROJECT_NAME}.x" PUBLIC OpenMP::OpenMP_CXX)
else()
    message("OpenMP not found, continuing without it...")
endif()

find_package(GSL REQUIRED)
target_link_libraries("${PROJECT_NAME}.x" PUBLIC GSL::gsl)

message("Build type       : " ${CMAKE_BUILD_TYPE}) #Prints Debug or Release, so it's clear if the code has actually been compiled with/without optimizations
message("Eigen dir        : " ${EIGEN3_INCLUDE_DIR}) #there actually is a more modern way using find_package but that did not work for me
message("Boost dir        : " ${BOOST_INCLUDE_DIR})