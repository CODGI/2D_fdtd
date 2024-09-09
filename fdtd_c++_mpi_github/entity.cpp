#include "entity.h"
#include "vector"
#include "cmath"
#include <iostream>


using namespace std;


int entity::getNx() {
    return Nx;
}

int entity::getNy() {
    return Ny;
}

rectangle::rectangle(double lx, double ly) {
    this->lx = lx;
    this->ly = ly;
}

vector<vector<bool> > rectangle::structure(double dx, double dy) {
    this->Nx = (int) (this->lx/dx);
    this->Ny = (int) (this->ly/dy);
    return vector<vector<bool> >(this->Nx, vector<bool>(this->Ny,true));
}




ellipse::ellipse(double rx, double ry) {
    this->rx = rx;
    this->ry = ry;
}

vector<vector<bool> > ellipse::structure(double dx, double dy) {
    this->Nx = (int) (2*(this->rx)/dx);
    this->Ny = (int) (2*(this->ry)/dy);
    if (this->Nx%2 == 0) { this->Nx = this->Nx+1; };
    if (this->Ny%2 == 0) { this->Ny = this->Ny+1; };
    vector<vector<bool> > occ = vector<vector<bool> >(this->Nx, vector<bool>(this->Ny,false));
    int midx = (int) ((this->Nx-1)/2);
    int midy = (int) ((this->Ny-1)/2);
    for (int i=0; i < occ.size(); ++i) {
        for (int j=0; j< occ[i].size();++j) {
            double x = ((i-midx)*dx);
            double y = ((j-midy)*dy);
            if (pow(x/this->rx,2)+pow(y/this->ry,2) <= 1) 
            {
                occ[i][j] = 1;
            }
        } 
    }
    return occ;
}