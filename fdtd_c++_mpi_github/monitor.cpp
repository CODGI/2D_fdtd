#include "monitor.h"
#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

monitor::monitor(double Lx, double Ly) {
    Lx = Lx;
    Ly = Ly;
    cerr << Lx << "," << Ly << endl;
}

void monitor::printToCSV(vector<vector<double>> v, int number, double t) {
    cerr << number << "," << t << endl;
    for (auto& row:v) {
        for (auto& col:row) {    
            cerr << col << ",";
        }
        cerr << endl;
    }    
}