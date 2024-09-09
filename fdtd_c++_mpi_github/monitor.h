#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
using namespace std;

class monitor {
    public:
        monitor(double, double);
        void printToCSV(vector<vector<double>>, int, double);
    private:
        double Lx;
        double Ly;
};