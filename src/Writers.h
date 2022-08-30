//
// Created by machiels on 30/08/2022.
//

#ifndef HOTPLASMAKERNEL_WRITERS_H
#define HOTPLASMAKERNEL_WRITERS_H

#include <string>
#include "Eigen/Dense"

class Writers {
public:
    Writers();

    /**
     * Write column data to csv file
     * @param fileName
     * @param dat
     */
    static void writeCsv(const std::string& fileName, Eigen::VectorXcd &dat);

    static void writeCsv(const std::string &fileName, std::complex<double> *dat, int size);

    static void writeCsv(const std::string &fileName, double *dat, int size);

    static void writeCsv2(const std::string& fileName, Eigen::MatrixXcd &dat, bool isReal);

};


#endif //HOTPLASMAKERNEL_WRITERS_H
