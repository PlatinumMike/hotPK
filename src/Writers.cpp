//
// Created by machiels on 30/08/2022.
//

#include "Writers.h"
#include <fstream>
#include <iomanip>

void Writers::writeCsv(const std::string& fileName, Eigen::VectorXcd &dat){
    //todo: write to HDF5
    std::ofstream csvFile;
    csvFile.open(fileName);
    for (Eigen::Index i=0; i<dat.size();i++){
        csvFile << std::setprecision (17) << dat(i).real() <<',' << dat(i).imag()<<std::endl;
    }
    csvFile.close();
}

void Writers::writeCsv(const std::string& fileName, std::complex<double> *dat, int size){
    //todo: write to HDF5
    std::ofstream csvFile;
    csvFile.open(fileName);
    for (Eigen::Index i=0; i<size;i++){
        csvFile << std::setprecision (17) << dat[i].real() <<',' << dat[i].imag()<<std::endl;
    }
    csvFile.close();
}

void Writers::writeCsv(const std::string& fileName, double *dat, int size){
    //todo: write to HDF5
    std::ofstream csvFile;
    csvFile.open(fileName);
    for (Eigen::Index i=0; i<size;i++){
        csvFile << std::setprecision (17) << dat[i] <<std::endl;
    }
    csvFile.close();
}

void Writers::writeCsv2(const std::string& fileName, Eigen::MatrixXcd &dat, bool isReal) {
    //todo: write to HDF5
    std::ofstream csvFile;
    csvFile.open(fileName);
    if (isReal) {
        //write real part
        for (Eigen::Index i = 0; i < dat.rows(); i++) {
            for (Eigen::Index j = 0; j < dat.cols() - 1; j++) {
                csvFile << std::setprecision(17) << dat(i, j).real() << ',';
            }
            csvFile << std::setprecision(17)
                    << dat(i, dat.cols() - 1).real(); //write separately to avoid trailing comma
            csvFile << std::endl;
        }
    } else {
        //write imag part
        for (Eigen::Index i = 0; i < dat.rows(); i++) {
            for (Eigen::Index j = 0; j < dat.cols() - 1; j++) {
                csvFile << std::setprecision(17) << dat(i, j).imag() << ',';
            }
            csvFile << std::setprecision(17) << dat(i, dat.cols() - 1).imag();
            csvFile << std::endl;
        }
    }
    csvFile.close();
}

Writers::Writers() = default;