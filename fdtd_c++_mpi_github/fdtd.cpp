#include "fdtd.h"
#include <vector>
#include "domainHandler.h"
#include "source.h"
#include <string.h>
#include <iostream>
using namespace std;


fdtd::fdtd(domainHandler& dh, double dt, double T, int sample,  source* J, monitor* m) {
    this->Nx = dh.Nx;
    this->Ny = dh.Ny;
    this->dx = dh.dx;
    this->dy = dh.dy;
    this->dt = dt;
    this->T = T;
    this->sample = sample;
    this->J = J;
    this->m = m;
    N_LorentzPoles = dh.N_LorentzPoles;
    N_DrudePoles = dh.N_DrudePoles;
    vector<vector<vector<double>>> eps_mu_staggered = vector<vector<vector<double>>>(8,vector<vector<double>>(Nx, vector<double>(Ny,0)));
    eps_mu_staggered[0] = dh.eps; eps_mu_staggered[1] = dh.eps; eps_mu_staggered[2] = dh.eps; eps_mu_staggered[3] = dh.eps; 
    eps_mu_staggered[4] = dh.mu; eps_mu_staggered[5] = dh.mu; eps_mu_staggered[6] = dh.mu; eps_mu_staggered[7] = dh.mu; 
    //eps_mu_staggered = interpolateToGrid(eps_mu_staggered);
    vector<vector<vector<double>>> sigmaStaggered = vector<vector<vector<double>>>(8,vector<vector<double>>(Nx, vector<double>(Ny,0)));
    sigmaStaggered[0] = dh.sigmaEy; sigmaStaggered[1] = dh.sigmaEx; sigmaStaggered[2] = dh.sigmaEx; sigmaStaggered[3] = dh.sigmaEy; 
    sigmaStaggered[4] = dh.sigmaHy; sigmaStaggered[5] = dh.sigmaHx; sigmaStaggered[6] = dh.sigmaHx; sigmaStaggered[7] = dh.sigmaHy; 
    //sigmaStaggered = interpolateToGrid(sigmaStaggered);
    generateDrudeLorentzGrid(dh);
    generateCoefficients(eps_mu_staggered, sigmaStaggered);
    J_drude = vector<vector<vector<vector<complex<double>>>>>(4,vector<vector<vector<complex<double>>>>(N_DrudePoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)))));
    J_drude_prev = vector<vector<vector<vector<complex<double>>>>>(4,vector<vector<vector<complex<double>>>>(N_DrudePoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)))));
    J_Lorentz = vector<vector<vector<vector<complex<double>>>>>(4,vector<vector<vector<complex<double>>>>(N_LorentzPoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)))));
    J_Lorentz_prev = vector<vector<vector<vector<complex<double>>>>>(4,vector<vector<vector<complex<double>>>>(N_LorentzPoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)))));
    this->Ex = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ey = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ezx = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ezy = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ex_prev = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ey_prev = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ezx_prev = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Ezy_prev = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Hx = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Hy = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Hzx = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    this->Hzy = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0))); 
}

void fdtd::generateDrudeLorentzGrid(domainHandler& dh) {
    cout << "Drude poles:" << N_DrudePoles << endl;
    N_DrudePoles = dh.N_DrudePoles;
    DrudeA = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_DrudePoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    DrudeA[0] = dh.DrudeA; DrudeA[1] = dh.DrudeA; DrudeA[2] = dh.DrudeA; DrudeA[3] = dh.DrudeA;
    //DrudeA = interpolateToGrid(DrudeA, N_DrudePoles);
    DrudeB = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_DrudePoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    DrudeB[0] = dh.DrudeB; DrudeB[1] = dh.DrudeB; DrudeB[2] = dh.DrudeB; DrudeB[3] = dh.DrudeB;
    //DrudeB = interpolateToGrid(DrudeB, N_DrudePoles);
    DrudeC = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_DrudePoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    DrudeC[0] = dh.DrudeC; DrudeC[1] = dh.DrudeC; DrudeC[2] = dh.DrudeC; DrudeC[3] = dh.DrudeC;
    //DrudeC = interpolateToGrid(DrudeC, N_DrudePoles);

    N_LorentzPoles = dh.N_LorentzPoles;
    LorentzA = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_LorentzPoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    LorentzA[0] = dh.LorentzA; LorentzA[1] = dh.LorentzA; LorentzA[2] = dh.LorentzA; LorentzA[3] = dh.LorentzA;
    //LorentzA = interpolateToGrid(LorentzA, N_LorentzPoles);
    LorentzB = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_LorentzPoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    LorentzB[0] = dh.LorentzB; LorentzB[1] = dh.LorentzB; LorentzB[2] = dh.LorentzB; LorentzB[3] = dh.LorentzB;
    //LorentzB = interpolateToGrid(LorentzB, N_LorentzPoles);
    LorentzC = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N_LorentzPoles, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    LorentzC[0] = dh.LorentzC; LorentzC[1] = dh.LorentzC; LorentzC[2] = dh.LorentzC; LorentzC[3] = dh.LorentzC;
    //LorentzC = interpolateToGrid(LorentzC, N_LorentzPoles);
}

void fdtd::generateCoefficients(vector<vector<vector<double>>>& eps_mu_staggered, vector<vector<vector<double>>>& sigma_staggered) {
    coefficient1 = vector<vector<vector<double>>>(8,vector<vector<double>>(Nx, vector<double>(Ny, 0)));
    coefficient2 = vector<vector<vector<double>>>(8,vector<vector<double>>(Nx, vector<double>(Ny, 0)));
    coefficient3 = vector<vector<vector<double>>>(8,vector<vector<double>>(Nx, vector<double>(Ny, 0)));
    for (int k=0; k < 8; k++) {
        for (int i=0; i < Nx; i++) {
            for (int j=0; j < Ny; j++) {
                double g = 0;
                if (k<4) {
                    for (int l=0; l<N_DrudePoles; l++) {
                        g += DrudeC[k][l][i][j]/4;
                    }
                    for (int l=0; l<N_LorentzPoles; l++) {
                        g += LorentzC[k][l][i][j]/4;
                    }
                }
                coefficient1[k][i][j] = dt/(eps_mu_staggered[k][i][j] + 0.5*dt*sigma_staggered[k][i][j]+g);
                coefficient2[k][i][j] = (eps_mu_staggered[k][i][j] - 0.5*dt*sigma_staggered[k][i][j])/(eps_mu_staggered[k][i][j] + 0.5*dt*sigma_staggered[k][i][j]+g);
                coefficient3[k][i][j] = g/(eps_mu_staggered[k][i][j] + 0.5*dt*sigma_staggered[k][i][j]+g);
            }
        }
    }
}

vector<vector<vector<vector<double>>>> fdtd::interpolateToGrid(vector<vector<vector<vector<double>>>>& a, int N) {
    vector<vector<vector<vector<double>>>> newA = vector<vector<vector<vector<double>>>>(4,vector<vector<vector<double>>>(N, vector<vector<double>>(Nx,vector<double>(Ny,0))));
    for (int k=0; k < N;k++) {
        for (int i=0; i < Nx-1; i++) {
            for (int j=0; j  < Ny-1; j++) {
                newA[0][k][i][j] = (a[0][k][i+1][j]+a[0][k][i][j])/2;
                newA[1][k][i][j] = (a[1][k][i][j+1]+a[1][k][i][j])/2;
                newA[2][k][i][j] = a[2][k][i][j];
                newA[3][k][i][j] = a[3][k][i][j];
            }
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::deriv_Yee(vector<vector<complex<double>>>& arr, string axis, string field, double dx) {
    vector<vector<complex<double>>> newArr = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    if (axis=="x") {
        if (field=="H") {
            for (int i=1; i<arr.size()-1; i++) {
                for (int j=1; j<arr[i].size()-1;j++) {
                    newArr[i][j] = (arr[i][j]-arr[i-1][j])/dx;
                }
            }
        } else if (field=="E") {
            for (int i=1; i<arr.size()-1; i++) {
                for (int j=1; j<arr[i].size()-1;j++) {
                    newArr[i][j] = (arr[i+1][j]-arr[i][j])/dx;
                }
            }
        } else { 
            cout << "Field has to be E or H" << endl;
            throw exception();
            }
    }
    if (axis=="y") {
        if (field=="H") {
            for (int i=1; i<arr.size()-1; i++) {
                for (int j=1; j<arr[i].size()-1;j++) {
                    newArr[i][j] = (arr[i][j]-arr[i][j-1])/dx;;
                }
            }
        } else if (field=="E") {
            for (int i=1; i<arr.size()-1; i++) {
                for (int j=1; j<arr[i].size()-1;j++) {
                    newArr[i][j] = (arr[i][j+1]-arr[i][j])/dx;
                }
            }
        } else { 
            cout << "Field has to be E or H" << endl;
            throw exception();
            }
    }
    return newArr;
}

vector<vector<complex<double>>> fdtd::generateDispCurrent(int index) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx, vector<complex<double>>(Ny,0));
    for (int i=0; i< Nx; i++) {
        for (int j=0; j<Ny; j++) {
            for (int l=0; l<N_DrudePoles; l++) {
                if (DrudeA[index][l][i][j] < 0) {
                cout << "A: " << DrudeA[index][l][i][j]<< endl;
                }
                if (DrudeB[index][l][i][j] < 0) {
                cout << "B: " << DrudeB[index][l][i][j]<< endl;
                }
                if (DrudeC[index][l][i][j] < 0) {
                cout << "C: " << DrudeC[index][l][i][j]<< endl;
                }
                newA[i][j] += 0.5*(1.0+DrudeA[index][l][i][j])*J_drude[index][l][i][j];
                newA[i][j] += 0.5*DrudeB[index][l][i][j]*J_drude_prev[index][l][i][j];
            }
            for (int l=0; l<N_LorentzPoles; l++) {
                newA[i][j] += 0.5*(1.0+LorentzA[index][l][i][j]*J_Lorentz[index][l][i][j]);
                newA[i][j] += 0.5*LorentzA[index][l][i][j]*J_Lorentz_prev[index][l][i][j];
            }
        }
    }
    return newA;
}

void fdtd::updateCurrent(vector<vector<complex<double>>>& E_new, vector<vector<complex<double>>>& E_prev, int index) {
    vector<vector<vector<complex<double>>>> J_Drude_new = vector<vector<vector<complex<double>>>>(N_DrudePoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0))));
    vector<vector<vector<complex<double>>>> J_Lorentz_new = vector<vector<vector<complex<double>>>>(N_LorentzPoles,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0))));
    for (int i=0; i < Nx; i++) {
        for (int j=0; j<Ny; j++) {
            complex<double> E_change = (E_new[i][j] - E_prev[i][j])/(2*dt);
            for (int l=0; l < N_DrudePoles; l++) {
                J_Drude_new[l][i][j] = DrudeA[index][l][i][j] * J_drude[index][l][i][j] + DrudeB[index][l][i][j] * J_drude_prev[index][l][i][j]+ DrudeC[index][l][i][j]*E_change;
            }
            for (int l=0; l < N_LorentzPoles; l++) {
                J_Lorentz_new[l][i][j] = LorentzA[index][l][i][j] * J_Lorentz[index][l][i][j] + LorentzB[index][l][i][j] * J_Lorentz_prev[index][l][i][j]+ LorentzC[index][l][i][j]*E_change;
            }
        }
    }
    J_drude_prev[index] = J_drude[index];
    J_drude[index] = J_Drude_new;
    J_Lorentz_prev[index] = J_Lorentz[index];
    J_Lorentz[index] = J_Lorentz_new;
}


void fdtd::makeStep(double t) {
    vector<vector<vector<complex<double>>>> update = vector<vector<vector<complex<double>>>>(8,vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0))));
    vector<vector<vector<complex<double>>>> j = J->radiate(t);
    
    //Ex
    vector<vector<complex<double>>> inEx = add(Hzx, Hzy);
    inEx = deriv_Yee(inEx,"y", "H", dy);
    inEx = subtract(inEx, j[0]);
    update[0] = multiply(inEx,coefficient1[0]);
    vector<vector<complex<double>>> inEx2 = multiply(Ex,coefficient2[0]);
    update[0] = add(update[0], inEx2);
    vector<vector<complex<double>>> inEx3 = multiply(Ex_prev,coefficient3[0]);
    update[0] = add(update[0], inEx3);
    vector<vector<complex<double>>> inEx4 = generateDispCurrent(0);
    inEx4 = multiply(inEx4, coefficient1[0]);
    update[0] = subtract(update[0], inEx4);
    updateCurrent(update[0], Ex_prev, 0);
    Ex_prev = Ex;
    Ex = update[0];

    //Ey
    vector<vector<complex<double>>> inEy = add(Hzx, Hzy);
    inEy = deriv_Yee(inEy,"x", "H", dx);
    inEy = scale(inEy,-1);
    inEy = subtract(inEy, j[1]);
    update[1] = multiply(inEy,coefficient1[1]);
    vector<vector<complex<double>>> inEy2 = multiply(Ey,coefficient2[1]);
    update[1] = add(update[1], inEy2);
    vector<vector<complex<double>>> inEy3 = multiply(Ey_prev,coefficient3[1]);
    update[1] = add(update[1], inEy3);
    vector<vector<complex<double>>> inEy4 = generateDispCurrent(1);
    inEy4 = multiply(inEy4, coefficient1[1]);
    update[1] = subtract(update[1], inEy4);
    updateCurrent(update[1], Ey_prev, 1);
    Ey_prev = Ey;
    Ey = update[1];

    //Ezx
    vector<vector<complex<double>>> inEzx = Hy;
    inEzx = deriv_Yee(inEzx,"x", "H", dx);
    inEzx = subtract(inEzx, j[2]);
    update[2] = multiply(inEzx,coefficient1[2]);
    vector<vector<complex<double>>> inEzx2 = multiply(Ezx,coefficient2[2]);
    update[2] = add(update[2], inEzx2);
    vector<vector<complex<double>>> inEzx3 = multiply(Ezx_prev,coefficient3[2]);
    update[2] = add(update[2], inEzx3);
    vector<vector<complex<double>>> inEzx4 = generateDispCurrent(2);
    inEzx4 = multiply(inEzx4, coefficient1[2]);
    update[2] = subtract(update[2], inEzx4);
    updateCurrent(update[2], Ezx_prev, 1);
    Ezx_prev = Ezx;
    Ezx = update[2];

    //Ezy
    vector<vector<complex<double>>> inEzy = Hx;
    inEzy = deriv_Yee(inEzy,"y", "H", dx);
    inEzy = scale(inEzy,-1);
    inEzy = subtract(inEzy, j[3]);
    update[3] = multiply(inEzy,coefficient1[3]);
    vector<vector<complex<double>>> inEzy2 = multiply(Ezy,coefficient2[3]);
    update[3] = add(update[3], inEzy2);
    vector<vector<complex<double>>> inEzy3 = multiply(Ezy_prev,coefficient3[3]);
    update[3] = add(update[3], inEzy3);
    vector<vector<complex<double>>> inEzy4 = generateDispCurrent(3);
    inEzy4 = multiply(inEzy4, coefficient1[3]);
    update[3] = subtract(update[3], inEzy4);
    updateCurrent(update[3], Ezy_prev, 1);
    Ezy_prev = Ezy;
    Ezy = update[3];

    //Hx
    vector<vector<complex<double>>> inHx = add(Ezx, Ezy);
    inHx = deriv_Yee(inHx,"y", "E", dy);
    inHx = scale(inHx, -1);
    update[4] = multiply(inHx, coefficient1[4]);
    vector<vector<complex<double>>> inHx2 = multiply(Hx, coefficient2[4]);
    update[4] = add(update[4], inHx2);
    Hx = update[4];

    //Hy
    vector<vector<complex<double>>> inHy = add(Ezx, Ezy);
    inHy = deriv_Yee(inHy,"x", "E", dx);
    update[5] = multiply(inHy, coefficient1[5]);
    vector<vector<complex<double>>> inHy2 = multiply(Hy, coefficient2[5]);
    update[5] = add(update[5], inHy2);
    Hy = update[5];

    //Hzx
    vector<vector<complex<double>>> inHzx = Ey;
    inHzx = deriv_Yee(inHzx,"x", "E", dx);
    inHzx = scale(inHzx,-1);
    update[6] = multiply(inHzx, coefficient1[6]);
    vector<vector<complex<double>>> inHzx2 = multiply(Hzx, coefficient2[6]);
    update[6] = add(update[6], inHzx2);
    Hzx = update[6];

    //Hzy
    vector<vector<complex<double>>> inHzy = Ex;
    inHzy = deriv_Yee(inHzy,"y", "E", dy);
    update[7] = multiply(inHzy, coefficient1[7]);
    vector<vector<complex<double>>> inHzy2 = multiply(Hzy, coefficient2[7]);
    update[7] = add(update[7], inHzy2);
    Hzy = update[7];
    
    if ((this->step)%(this->sample) == 0) {
        //cout << t << endl;
        vector<vector<double>> I = calculateIntensity();
        m->printToCSV(I, step,t);
    }
}

void fdtd::run() {
    double t = 0;
    step = 0;
    while (t<this->T) {
        makeStep(t);
        cout << t << endl;
        t += dt;
        step++;
    }
}

vector<vector<vector<double>>> fdtd::interpolateToGrid(vector<vector<vector<double>>>& a) {
    vector<vector<vector<double>>> newA = vector<vector<vector<double>>>(a.size(),vector<vector<double>>(Nx, vector<double>(Ny,1)));
    for (int i=0; i<Nx-1; i++) {
        for (int j=0; j<Ny-1; j++) {
            newA[0][i][j] = (a[0][i+1][j]+a[0][i][j])/2;
            newA[1][i][j] = (a[1][i][j+1]+a[1][i][j])/2;
            newA[2][i][j] = a[2][i][j];
            newA[3][i][j] = a[3][i][j];
            newA[4][i][j] = (a[4][i][j+1]+a[4][i][j])/2;
            newA[5][i][j] = (a[5][i+1][j]+a[5][i][j])/2;
            newA[6][i][j] = (a[6][i][j]+a[6][i+1][j]+a[6][i][j+1]+a[6][i+1][j+1])/4;
            newA[7][i][j] = (a[7][i][j]+a[7][i+1][j]+a[7][i][j+1]+a[7][i+1][j+1])/4;
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::add(vector<vector<complex<double>>>& a1, vector<vector<complex<double>>>& a2) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = a1[i][j] + a2[i][j];
        }
    }
    return newA;
}

vector<vector<double>> fdtd::add(vector<vector<double>>& a1, vector<vector<double>>& a2) {
    vector<vector<double>> newA = vector<vector<double>>(Nx,vector<double>(Ny, double(0)));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = a1[i][j] + a2[i][j];
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::subtract(vector<vector<complex<double>>>& a1, vector<vector<complex<double>>>& a2) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = a1[i][j] - a2[i][j];
        }
    }
    return newA;
}

vector<vector<double>> fdtd::subtract(vector<vector<double>>& a1, vector<vector<double>>& a2) {
    vector<vector<double>> newA = vector<vector<double>>(Nx,vector<double>(Ny,0));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = a1[i][j] - a2[i][j];
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::multiply(vector<vector<complex<double>>>& a1, vector<vector<double>>& a2) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = complex<double>((a1[i][j]).real() * a2[i][j], (a1[i][j]).imag() * a2[i][j]);
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::divide(vector<vector<complex<double>>>& a1, vector<vector<double>>& a2) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = complex<double>((a1[i][j]).real() / a2[i][j], (a1[i][j]).imag() / a2[i][j]);
        }
    }
    return newA;
}

vector<vector<double>> fdtd::divide(vector<vector<double>>& a1, vector<vector<double>>& a2) {
    vector<vector<double>> newA = vector<vector<double>>(Nx,vector<double>(Ny, 0));
    for (int i=0; i<a1.size(); i++) {
        for (int j=0; j<a1[i].size();j++) {
            newA[i][j] = a1[i][j]/ a2[i][j];
        }
    }
    return newA;
}

vector<vector<complex<double>>> fdtd::scale(vector<vector<complex<double>>>& a, double DT) {
    vector<vector<complex<double>>> newA = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=0; i<a.size(); i++) {
        for (int j=0; j<a[i].size();j++) {
            newA[i][j] = complex<double>((a[i][j]).real() *DT, (a[i][j]).imag() * DT);
        }
    }
    return newA;
}

vector<vector<double>> fdtd::scale(vector<vector<double>>& a, double DT) {
    vector<vector<double>> newA = vector<vector<double>>(Nx,vector<double>(Ny, 0));
    for (int i=0; i<a.size(); i++) {
        for (int j=0; j<a[i].size();j++) {
            newA[i][j] = a[i][j]*DT;
        }
    }
    return newA;
}

vector<vector<double>> fdtd::calculateIntensity(){
    vector<vector<complex<double>>> Ex_interpolated = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    vector<vector<complex<double>>> Ey_interpolated = vector<vector<complex<double>>>(Nx,vector<complex<double>>(Ny, complex<double>(0,0)));
    for (int i=1; i<Nx;i++) {
        for (int j=1; j<Ny;j++){
            Ex_interpolated[i][j] = (Ex[i][j]+Ex[i-1][j])/2.0;
            Ey_interpolated[i][j] = (Ey[i][j]+Ey[i][j-1])/2.0;
        }
    }
    vector<vector<double>> Ex2 = squareMatrix(Ex_interpolated);
    vector<vector<double>> Ey2 = squareMatrix(Ey_interpolated);
    vector<vector<complex<double>>> Ez = add(Ezx, Ezy);
    vector<vector<double>> Ez2 = squareMatrix(Ez);
    vector<vector<double>> I = add(Ex2, Ey2);
    I = add(I, Ez2);
    return I;
}

vector<vector<double>> fdtd::squareMatrix(vector<vector<complex<double>>>& v) {
    vector<vector<double>> vSquared = vector<vector<double>>(Nx,vector<double>(Ny, 0));
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j<Ny; j++) {
            vSquared[i][j] = pow(std::abs(v[i][j]),2);
        }
    }
    return vSquared;
}