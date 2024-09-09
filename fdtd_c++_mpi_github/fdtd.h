#pragma once

#include <vector>
#include "domainHandler.h"
#include "source.h"
#include "string.h"
#include "monitor.h"

class fdtd {
    public:
        fdtd(domainHandler&, double, double, int, source*, monitor*);
        void run();
    private:
        int Nx;
        int Ny;
        double dx;
        double dy;
        double dt;
        double T;
        int sample;
        int step;
        int fourierPoints;
        double fourierWidth;
        double fourierCenter;
        vector<double> wavelengths; 
        vector<vector<double>> sigmaEx;
        vector<vector<double>> sigmaEy;
        vector<vector<double>> sigmaHx;
        vector<vector<double>> sigmaHy;
        source* J;
        monitor* m;
        vector<vector<complex<double>>> Ex;
        vector<vector<complex<double>>> Ey;
        vector<vector<complex<double>>> Ezx;
        vector<vector<complex<double>>> Ezy;
        vector<vector<complex<double>>> Ex_prev;
        vector<vector<complex<double>>> Ey_prev;
        vector<vector<complex<double>>> Ezx_prev;
        vector<vector<complex<double>>> Ezy_prev;
        vector<vector<vector<vector<complex<double>>>>> J_drude;
        vector<vector<vector<vector<complex<double>>>>> J_drude_prev;
        vector<vector<vector<vector<complex<double>>>>> J_Lorentz;
        vector<vector<vector<vector<complex<double>>>>> J_Lorentz_prev;
        vector<vector<vector<double>>> coefficient1;
        vector<vector<vector<double>>> coefficient2;
        vector<vector<vector<double>>> coefficient3;
        vector<vector<complex<double>>> Hx;
        vector<vector<complex<double>>> Hy;
        vector<vector<complex<double>>> Hzx;
        vector<vector<complex<double>>> Hzy;

        int N_LorentzPoles;
        int N_DrudePoles;
        void generateDrudeLorentzGrid(domainHandler&);
        void generateCoefficients(vector<vector<vector<double>>>&, vector<vector<vector<double>>>&);
        vector<vector<vector<vector<double>>>> DrudeA;
        vector<vector<vector<vector<double>>>> DrudeB;
        vector<vector<vector<vector<double>>>> DrudeC;
        vector<vector<vector<vector<double>>>> LorentzA;
        vector<vector<vector<vector<double>>>> LorentzB;
        vector<vector<vector<vector<double>>>> LorentzC;

        vector<vector<vector<double>>> interpolateToGrid(vector<vector<vector<double>>>&);
        vector<vector<vector<vector<double>>>> interpolateToGrid(vector<vector<vector<vector<double>>>>&, int);

        vector<vector<complex<double>>> deriv_Yee(vector<vector<complex<double>>>&, string, string, double);
        void makeStep(double);
        vector<vector<complex<double>>> generateDispCurrent(int);
        void updateCurrent(vector<vector<complex<double>>>&, vector<vector<complex<double>>>& , int);

        vector<vector<complex<double>>> add(vector<vector<complex<double>>>&, vector<vector<complex<double>>>&);
        vector<vector<double>> add(vector<vector<double>>& a1, vector<vector<double>>& a2);
        vector<vector<complex<double>>> subtract(vector<vector<complex<double>>>&, vector<vector<complex<double>>>&);
        vector<vector<double>> subtract(vector<vector<double>>&, vector<vector<double>>&);
        vector<vector<complex<double>>> multiply(vector<vector<complex<double>>>&, vector<vector<double>>&);
        vector<vector<complex<double>>> divide(vector<vector<complex<double>>>&, vector<vector<double>>&);
        vector<vector<double>> divide(vector<vector<double>>&, vector<vector<double>>&);
        vector<vector<complex<double>>> scale(vector<vector<complex<double>>>&, double);
        vector<vector<complex<double>>> scale(vector<vector<complex<double>>>&, complex<double>);
        vector<vector<double>> scale(vector<vector<double>>&, double);

        vector<vector<double>> calculateIntensity();
        vector<vector<double>> squareMatrix(vector<vector<complex<double>>>& v);


};