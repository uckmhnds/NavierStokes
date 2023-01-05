#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED

#include <cmath>
extern int nodeNumber, cellNumber;

int nodeNumber  = 948;
int cellNumber  = 1800;


extern char filename[];

char filename[] = "naca.dat";


extern int timestep, maxstep;

int timestep    = 1;
int maxstep     = 1000;

extern float gamma_, gamm1_, PI, rk3cf[3];

float gamma_    = 1.4;
float gamm1_    = 0.4;
float PI        = 3.14;
float rk3cf[3]  = {0.3333, 0.5, 1.};

extern float fsmach, alpha, alphar, cfl;

float fsmach    = 0.3;
float alpha     = 0.;
float alphar    = alpha * PI / 180.;
float cfl       = 0.5;

extern float rho_, u_, v_, T_, p_, e_;

float rho_      = 1.;
float u_        = fsmach * cos(alphar);
float v_        = fsmach * sin(alphar);
float T_        = 1.;
float p_        = rho_ * T_ / gamma_;
float e_        = p_ / gamm1_ + 0.5 * rho_ * pow(fsmach, 2.);


float * dtmin   = new float[cellNumber];

#endif // VARIABLES_H_INCLUDED
