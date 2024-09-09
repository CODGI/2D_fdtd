#include "domainHandler.h"
#include "entity.h"
#include "material.h"
#include "position.h"
#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;

void domainHandler::printDomain(){
    cout << "Printing to files!" << endl;

    cout << "Print eps" << endl;
    ofstream epsFile("eps.csv");
    epsFile << Lx << "," << Ly << "," << Nx << "," << Ny << endl;
    for (auto& row:eps) {
        for (auto& col:row) {    
            epsFile << col << ",";
        }
        epsFile << endl;
    }

    cout << "Print sigmaE" << endl;
    ofstream sigmaFile("sigmaE.csv");
    sigmaFile << Lx << "," << Ly << "," << Nx << "," << Ny << endl;
    for (int i=0; i<sigmaEx.size(); i++) {
        for (int j=0; j<sigmaEx[i].size(); j++) {    
            sigmaFile << (sigmaEx[i][j]+sigmaEy[i][j]) << ",";
        }
        sigmaFile << endl;
    }

    cout << "Print mu" << endl;
    ofstream muFile("mu.csv");
    muFile << Lx << "," << Ly << "," << Nx << "," << Ny << endl;
    for (auto& row:mu) {
        for (auto& col:row) {    
            muFile << col << ",";
        }
        muFile << endl;
    }

    cout << "Print sigmaH" << endl;
    ofstream sigmaHFile("sigmaH.csv");
    sigmaHFile << Lx << "," << Ly << "," << Nx << "," << Ny << endl;
    for (int i=0; i<sigmaHx.size(); i++) {
        for (int j=0; j<sigmaHx[i].size(); j++) {    
            sigmaHFile << (sigmaHx[i][j]+sigmaHy[i][j]) << ",";
        }
        sigmaHFile << endl;
    }

}

domainHandler::domainHandler(double lx, double ly, double dx_t, double dy_t, double DT) {
    Lx = lx;
    Ly = ly;
    dt = DT;
        
    Nx = (int) (lx/dx_t);
    Ny = (int) (ly/dy_t);

    if (Nx%2 == 0) {
        Nx++;
    }
    if (Ny%2 == 0) {
        Ny++;
    }

    dx = (lx/Nx);
    dy = (ly/Ny);

    occ = vector<vector<bool> >(Nx, vector<bool>(Ny,false));
    eps = vector<vector<double> >(Nx, vector<double>(Ny,1));
    mu = vector<vector<double> > (Nx, vector<double>(Ny,1));
    sigmaEx = vector<vector<double > > (Nx, vector<double >(Ny,0));
    sigmaEy = vector<vector<double > > (Nx, vector<double >(Ny,0));
    sigmaHx = vector<vector<double > > (Nx, vector<double >(Ny,0));
    sigmaHy = vector<vector<double > > (Nx, vector<double >(Ny,0));

    N_LorentzPoles = 0;
    N_DrudePoles = 0;

}

void domainHandler::addPML(double sigma, double p) {
    int pmlWidthInCellsX = (int) ((p/2)*Nx);
    int pmlWidthInCellsY = (int) ((p/2)*Ny);
    for (int i = 0; i<pmlWidthInCellsX; i++) {
        for (int j=0; j<Ny;j++) {
            sigmaEx[i][j] = eps[i][j]*sigma*(((double)(pmlWidthInCellsX-i))/pmlWidthInCellsX);
            sigmaEx[Nx-1-i][j] = eps[Nx-1-i][j]*sigma*((double)(pmlWidthInCellsX-i)/pmlWidthInCellsX);
            sigmaHx[i][j] = sigmaEx[i][j]/eps[i][j];
            sigmaHx[Nx-1-i][j] = sigmaEx[Nx-1-i][j]/eps[Nx-1-i][j];
        }
        for (int j=0; j<Nx;j++) {
            sigmaEy[j][i] = eps[j][i]*sigma*(((double)(pmlWidthInCellsY-i))/pmlWidthInCellsY);
            sigmaEy[j][Ny-1-i] = eps[j][Ny-1-i]*sigma*((double)(pmlWidthInCellsY-i)/pmlWidthInCellsY);
            sigmaHy[j][i] = sigmaEy[j][i]/eps[j][i];
            sigmaHy[j][Ny-1-i] = sigmaEy[j][Ny-1-i]/eps[j][Ny-1-i];
        }
    }
}

void domainHandler::addMaterial(vector<vector<bool> > occ, material m) {
    for (int i=0; i<occ.size(); i++) {
        for(int j=0; j<occ[i].size(); j++){
            if (occ[i][j]) {
                this->eps[i][j] = m.eps;
                this->mu[i][j] = m.mu;
                this->sigmaEx[i][j] = m.sigmaE;
                this->sigmaEy[i][j] = m.sigmaE;
                this->sigmaHx[i][j] = m.sigmaH;
                this->sigmaHy[i][j] = m.sigmaH;

                if (m.N_LorentzPoles != 0) {
                    if (N_LorentzPoles < m.N_LorentzPoles) {
                        for (int n=N_LorentzPoles; n<m.N_LorentzPoles; n++) {
                            LorentzA[n] = vector<vector<double>>(Nx, vector<double>(Ny, 0));
                            LorentzB[n] = vector<vector<double>>(Nx, vector<double>(Ny, 0));
                            LorentzC[n] = vector<vector<double>>(Nx, vector<double>(Ny, 0));
                        }
                        N_LorentzPoles = m.N_LorentzPoles;
                    }
                    for (int k=0; k<N_LorentzPoles; k++) {
                        LorentzA[k][i][j] = (2.0-(m.LorentzOmegaP[k]*m.LorentzOmegaP[k]*dt*dt))/(m.LorentzDeltaP[k]*dt+1.0);
                        LorentzB[k][i][j] = (m.LorentzDeltaP[k]*dt - 1.0)/(m.LorentzDeltaP[k]*dt + 1.0);
                        LorentzC[k][i][j] = (m.LorentzDeltaEpsP[k]*m.LorentzOmegaP[k]*m.LorentzOmegaP[k]*dt*dt)/(m.LorentzDeltaP[k]*dt + 1.0);
                    }
                }

                if (m.N_DrudePoles != 0) {
                    if (N_DrudePoles < m.N_DrudePoles) {
                        for (int n=N_DrudePoles; n<m.N_DrudePoles; n++) {
                            DrudeA.push_back(vector<vector<double>>(Nx, vector<double>(Ny, 0)));
                            DrudeB.push_back(vector<vector<double>>(Nx, vector<double>(Ny, 0)));
                            DrudeC.push_back(vector<vector<double>>(Nx, vector<double>(Ny, 0)));
                        }
                        N_DrudePoles = m.N_DrudePoles;
                        cout << "N_DrudePoles (domainHandler)" << N_DrudePoles << endl;
                    }
                    for (int k=0; k<N_DrudePoles; k++) {
                        DrudeA[k][i][j] = 2.0/(1.0 + 0.5*m.DrudeGammaQ[k]*dt);
                        DrudeB[k][i][j] = -(1.0 - 0.5*m.DrudeGammaQ[k]*dt)/(1.0 + 0.5*m.DrudeGammaQ[k]*dt);
                        DrudeC[k][i][j] = (m.DrudeOmegaQ[k]*m.DrudeOmegaQ[k]*dt*dt)/(1.0 + 0.5*m.DrudeGammaQ[k]*dt);
                    }
                }
            }

        }
    }
}

void domainHandler::addEntity(entity& e, position p, material m, bool debug = false, bool override = true) {
    vector< vector<bool> >occ_loc = e.structure(this->dx, this->dy);
    int Nx_loc = e.getNx();
    int Ny_loc = e.getNy();
    int pos_x_index = (int) (((this->Nx)-1)/2)+ (int) (p.x/(this->dx));
    int pos_y_index = (int) (((this->Ny)-1)/2)+ (int) (p.y/(this->dy));
    int smallest_x_index = pos_x_index-(int) Nx_loc/2;
    int smallest_y_index = pos_y_index-(int) Ny_loc/2;
    if (debug) {
        cout << "pos_x_index: " << pos_x_index << endl;
        cout << "pos_y_index: " << pos_y_index << endl;
        cout << "smallest_x_index: " << smallest_x_index << endl;
        cout << "smallest_y_index: " << smallest_y_index << endl;
        cout << "Nx_loc: " << Nx_loc << endl;
        cout << "Ny_loc: " << Ny_loc << endl;
    }
    if ((smallest_x_index < 0) || (smallest_y_index < 0)) {
        cout << "Object out of bound!" << endl;
        throw exception();
    }
    if ((smallest_x_index + Nx_loc > this->Nx) || (smallest_y_index + Ny_loc > this->Ny)) {
        cout << "Object out of bound!" << endl;
        throw std::exception();
    }
    vector< vector<bool> >occ_ext = vector<vector<bool> >(this->Nx, vector<bool>(this->Ny,false));
    for (int i=0; i < occ_loc.size(); ++i) {
        for (int j=0; j< occ_loc[i].size();++j) {
            if (override) {
                occ_ext[smallest_x_index+i][smallest_y_index+j] = occ_loc[i][j];
                this->occ[smallest_x_index+i][smallest_y_index+j] = occ_loc[i][j];
            }
            else {
                if ((occ_loc[i][j]) && (this->occ[smallest_x_index+i][smallest_y_index+j])) {
                    cout << "Object overlap detected" << endl;
                    throw std::exception();
                }
                else {
                    occ_ext[smallest_x_index+i][smallest_y_index+j] = occ_loc[i][j];
                    this->occ[smallest_x_index+i][smallest_y_index+j] = occ_loc[i][j];
                }
            }
        }
    }
    this->addMaterial(occ_ext, m);
}
