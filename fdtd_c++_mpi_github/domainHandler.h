#pragma once

#include "vector"
#include "entity.h"
#include "material.h"
#include "position.h"
using namespace std;

class domainHandler
    {
        public:
            void printDomain();
            domainHandler(double, double, double, double, double);
            void addEntity(entity&, position, material, bool, bool);
            void addPML(double, double);
            int Nx;
            int Ny;
            double Lx;
            double Ly;
            double dx;
            double dy;
            double dt;
            vector<vector<bool>> occ;
            vector<vector<double>> eps;
            vector<vector<double>> mu;
            vector<vector<double>> sigmaEx;
            vector<vector<double>> sigmaEy;
            vector<vector<double>> sigmaHx;
            vector<vector<double>> sigmaHy;
            int N_LorentzPoles;
            int N_DrudePoles;
            vector<vector<vector<double>>> DrudeA;
            vector<vector<vector<double>>> DrudeB;
            vector<vector<vector<double>>> DrudeC;
            vector<vector<vector<double>>> LorentzA;
            vector<vector<vector<double>>> LorentzB;
            vector<vector<vector<double>>> LorentzC;
            void addMaterial(vector<vector<bool> >, material);
    };