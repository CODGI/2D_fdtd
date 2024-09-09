#pragma once

#include "vector"
using namespace std;

class entity 
{
    public:
        virtual vector< vector<bool> > structure(double, double) = 0;
        int getNx();
        int getNy();
    protected:
        int Nx;
        int Ny;
};

class rectangle : public entity 
{
    public:
        rectangle(double, double);
        vector< vector<bool> > structure(double, double);
    private:
        double lx;
        double ly;
};

class ellipse : public entity 
{
    public:
        ellipse(double, double);
        vector< vector<bool> > structure(double, double);
    private:
        double rx;
        double ry;
};
