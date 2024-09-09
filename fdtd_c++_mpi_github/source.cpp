#include "source.h"
#include "position.h"
#include <cmath>
#include <bits/stdc++.h>
using namespace std;

source::source(double Lx, double Ly, double dx, double dy) {
    this->Lx = Lx;
    this->Ly = Ly;
    this->Nx = (int) (Lx/dx);
    this->Ny = (int) (Ly/dy);
    if ((this->Nx%2) == 0) {
        (this->Nx)++;
    }
    if ((this->Ny%2) == 0) {
        (this->Ny)++;
    }
    this->dx = (this->Lx)/(this->Nx);
    this->dy = (this->Ly)/(this->Ny);
}

dipoleSource::dipoleSource(double Lx, double Ly, double dx, double dy, 
                            position pos, double omega, vector<int> orientation, 
                            double amplitude, double T) : source(Lx, Ly, dx, dy) {
    
    pos_x_index = (int)(((Nx)-1)/2)+(int)((pos.x)/(dx));
    pos_y_index = (int)(((Ny)-1)/2)+(int)((pos.y)/(dy));
    this->omega = omega;
    this->orientation = orientation;
    this->amplitude = amplitude;
    this->T = T;
};

vector<vector<vector<complex<double>>>> dipoleSource::radiate(double t) {
    vector<vector<vector<complex<double>>>> J = vector<vector<vector<complex<double>>>>(4, vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny,complex<double>(0,0))));
    for (int k=0; k<J.size(); k++) {
        if (orientation[k] == 1) {
                complex<double> z(exp(-t/T)*cos(omega*t), exp(-t/T)*sin(omega*t));
                //complex<double> z(cos(omega*t), sin(omega*t));
                J[k][pos_x_index][pos_y_index] = amplitude*z;
            }
        else {
            continue;
        }
    }
    return J;

}

lineSource::lineSource(double Lx, double Ly, double dx, double dy, 
                            position pos, double length, string o, double omega, vector<int> orientation, 
                            double amplitude, double T) : source(Lx, Ly, dx, dy) {
    pos_x_index = (int)(((Nx)-1)/2)+(int)((pos.x)/(dx));
    pos_y_index = (int)(((Ny)-1)/2)+(int)((pos.y)/(dy));
    this->omega = omega;
    this->orientation = orientation;
    this->amplitude = amplitude;
    this->T = T;
    this->o = o;
    if (o == "x") {
        int halfLengthCells = (int) ((length/dx)/2);
        leftIndex = pos_x_index - halfLengthCells;
        rightIndex = pos_x_index + halfLengthCells;
    } else if (o == "y") {
        int halfLengthCells = (int) (length/dy);
        leftIndex = pos_y_index - halfLengthCells;
        rightIndex = pos_y_index + halfLengthCells;
    }
};

vector<vector<vector<complex<double>>>> lineSource::radiate(double t) {
    vector<vector<vector<complex<double>>>> J = vector<vector<vector<complex<double>>>>(4, vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny,complex<double>(0,0))));
    for (int k=0; k<J.size(); k++) {
        if (orientation[k] == 1) {
                complex<double> z(exp(-t/T)*cos(omega*t), exp(-t/T)*sin(omega*t));
                if (o == "x") {
                    for (int i = leftIndex; i < rightIndex; i++) {
                        J[k][i][pos_y_index] = amplitude*z;  
                    }
                } else if (o == "y") {
                    for (int i=leftIndex; i<rightIndex; i++) {
                        J[k][pos_x_index][i] = amplitude*z;
                    }
                }
            }
        else {
            continue;
        }
    }
    return J;

}

gaussianSource::gaussianSource(double Lx, double Ly, double dx, double dy, 
                            position pos, double length, string o, double omega, vector<int> orientation, 
                            double amplitude, double T, double waist, double distance) : source(Lx, Ly, dx, dy) {
    pos_x_index = (int)(((Nx)-1)/2)+(int)((pos.x)/(dx));
    pos_y_index = (int)(((Ny)-1)/2)+(int)((pos.y)/(dy));
    this->omega = omega;
    this->orientation = orientation;
    this->T = T;
    this->o = o;
    double zR = 0.5*omega*pow(waist,2);
    this->wz = waist*sqrt(1+pow(distance/zR,2));
    this->Rz = distance*(1+pow(zR/distance,2));
    double phiz = atan(distance/zR);
    complex<double> phase(cos(phiz-omega*distance), sin(phiz-omega*distance));
    this->amplitude = amplitude*(waist/wz)*phase;
    if (o == "x") {
        int halfLengthCells = (int) ((length/dx)/2);
        leftIndex = pos_x_index - halfLengthCells;
        rightIndex = pos_x_index + halfLengthCells;
    } else if (o == "y") {
        int halfLengthCells = (int) (length/dy);
        leftIndex = pos_y_index - halfLengthCells;
        rightIndex = pos_y_index + halfLengthCells;
    }
};

vector<vector<vector<complex<double>>>> gaussianSource::radiate(double t) {
    vector<vector<vector<complex<double>>>> J = vector<vector<vector<complex<double>>>>(4, vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny,complex<double>(0,0))));
    for (int k=0; k<J.size(); k++) {
        if (orientation[k] == 1) {
                complex<double> z(exp(-t/T)*cos(omega*t), exp(-t/T)*sin(omega*t));
                if (o == "x") {
                    for (int i = leftIndex; i < rightIndex; i++) {
                        double r = dx*(i-pos_x_index);
                        double amp = exp(-pow(r/wz,2));
                        complex<double> phase_r(cos(-0.5*omega*pow(r,2)/Rz));
                        J[k][i][pos_y_index] = amplitude*amp*phase_r*z;  
                    }
                } else if (o == "y") {
                    for (int i=leftIndex; i<rightIndex; i++) {
                        double r = dx*(i-pos_y_index);
                        double amp = exp(-pow(r/wz,2));
                        complex<double> phase_r(cos(-0.5*omega*pow(r,2)/Rz));
                        J[k][i][pos_x_index] = amplitude*amp*phase_r*z;  
                    }
                }
            }
        else {
            continue;
        }
    }
    return J;

}