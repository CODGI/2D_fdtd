#pragma once

#include <complex>
#include <vector>
#include <string>
#include "position.h"
using namespace std;

class source
{
public:
    source(double, double, double, double);
    virtual vector<vector<vector<complex<double>>>> radiate(double) = 0;
    int Nx;
    int Ny;
    double Lx;
    double Ly;
    double dx;
    double dy;
};

class dipoleSource : public source
{
public:
    dipoleSource(double, double, double, double, position, double, vector<int>, double, double);
    vector<vector<vector<complex<double>>>> radiate(double);

private:
    int pos_x_index;
    int pos_y_index;
    double omega;
    vector<int> orientation;
    double amplitude;
    double T;
};

class lineSource : public source
{
public:
    lineSource(double, double, double, double, position, double, string, double, vector<int>, double, double);
    vector<vector<vector<complex<double>>>> radiate(double);

private:
    int leftIndex;
    int rightIndex;
    string o;
    int pos_x_index;
    int pos_y_index;
    double omega;
    vector<int> orientation;
    double amplitude;
    double T;
};

class gaussianSource : public source
{
    public:
        gaussianSource(double, double, double, double, position, double, string, double, vector<int>, double, double, double, double);
        vector<vector<vector<complex<double>>>> radiate(double);

    private: 
        int leftIndex;
        int rightIndex;
        string o;
        int pos_x_index;
        int pos_y_index;
        double omega;
        vector<int> orientation;
        complex<double> amplitude;
        double T;
        double Rz;
        double wz;
};