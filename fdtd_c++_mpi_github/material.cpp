#include "material.h"
#include "string"
#include <vector>
#include <iostream>
using namespace std;

material::material(string name, double eps, double mu, double sigmaE, double sigmaH) {
    this->name = name;
    this->eps = eps;
    this->mu = mu;
    this->sigmaE = sigmaE;
    this->sigmaH = sigmaH;
    N_DrudePoles = 0;
    N_LorentzPoles = 0;
}
void material::adddDrudePole(double omegaQ, double gammaQ) {
    DrudeOmegaQ.push_back(omegaQ);
    DrudeGammaQ.push_back(gammaQ);
    N_DrudePoles += 1;
    cout << "Drude poles (material) :"  << N_DrudePoles << endl;
}

void material::addLorentzPole(double deltaEps, double omegaP, double deltaP) {
    LorentzDeltaEpsP.push_back(deltaEps);
    LorentzOmegaP.push_back(omegaP);
    LorentzDeltaP.push_back(deltaP);
    N_LorentzPoles += 1;
}