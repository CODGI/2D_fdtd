#pragma once

#include "string"
#include <vector>
using namespace std;

class material {
    public:
        material(string, double, double, double,double);
        string name;
        double eps;
        double mu;
        double sigmaE;
        double sigmaH;
        int N_DrudePoles;
        int N_LorentzPoles;
        vector<double> LorentzOmegaP;
        vector<double> LorentzDeltaP;
        vector<double> LorentzDeltaEpsP;
        vector<double> DrudeOmegaQ;
        vector<double> DrudeGammaQ;
        void adddDrudePole(double, double);
        void addLorentzPole(double, double, double);
};
