#include <stdio.h>
#include <stdlib.h>
#include "mex.h"    
#include "matrix.h"   
#include <string.h> 
#include <math.h>
#include <time.h>
#include <iostream>

//------------- parameters --------------------------
//double C=1000, Iapp=0;
double C=1, Iapp;
double gL, gbarK, gbarCa, gbarh, gmi; // maximal conductances
//double EL=-40, EK, ECa=100, Eh=-20, Emi=-10; // reversal potentials
double EL=-50, EK=-80, ECa=100, Eh=-20, Emi=50; // reversal potentials (Skinner)
//double v1=0, v2=20; // for Ca current (Gabrielle)
double v1=0, v2=15; // for Ca current (Skinner)
//double v3=0, v4=15, phi=0.002; //for K current (Gabrielle)
double v3=0, v4=15, phi=0.0005; //for K current
//double v3=0, v4=15, phi=0.000002; //for K current
//double v5=78.3, v6=10.5, v7=-42.6, v8=87.3; // for H current
double v5=50, v6=10.5, v7=-42.6, v8=87.3; // for H current
//double vhalf_mi=-20, k_mi=-10;
double vhalf_mi=0, k_mi=-10;
//double V_th, V_slope=10.0, tau_s=10.0; // for Isyn
//double V_half=-20, V_slope=-2, tau_s=100, Esyn=-80, gsyn=10; // for Isyn
//double V_half=0, V_slope=0.001, tau_s=1, Esyn=-80, gsyn=6; // for Isyn
double V_half=0, V_slope=2, tau_s=1, Esyn=-80, gsyn=0.0; // for Isyn
double dt=0.02, TT, tgl; //step, total time, global time
int N, Neq=8; // number of steps, number of equations
double inp; // input
double *y, *f;
//---------------------------------------------------------
void system (void)
{    
// ML neuron #1
// Ca current
double minf = 0.5*(1+tanh((y[0]-v1)/v2));
//K current
double winf = 0.5*(1+tanh((y[0]-v3)/v4));
double lambdaw = phi*cosh((y[0]-v3)/(2*v4));
//H current
double h_inf = 1/(1+exp((y[0]+v5)/v6));
double tauh = 272-(-1499/(1+exp(-(y[0]+v7)/v8)));
// Imi current
double miinf=1/(1+exp((y[0]-vhalf_mi)/k_mi));

// Isyn
//double s_inf1, s_inf2;
//    if (y[0] > V_th) s_inf1= tanh( (y[0] - V_th) / V_slope);
//    else  s_inf1= 0;                                        
//double Isyn2 = gsyn*y[6]*(Esyn - y[3]);
    
//    if (y[3] > V_th) s_inf2= tanh( (y[3] - V_th) / V_slope);
//    else  s_inf2= 0;                                        
//double Isyn1 = gsyn*y[7]*(Esyn - y[0]); 

//double s_inf1 = 1.0/(1.0+exp((y[0]-V_half)/V_slope));
//double s_inf2 = 1.0/(1.0+exp((y[3]-V_half)/V_slope));	

// ML neuron #2
// Ca current
double minf2 = 0.5*(1+tanh((y[3]-v1)/v2));
//K current
double winf2 = 0.5*(1+tanh((y[3]-v3)/v4));
double lambdaw2 = phi*cosh((y[3]-v3)/(2*v4));
//H current
double h_inf2 = 1/(1+exp((y[3]+v5)/v6));
double tauh2 = 272-(-1499/(1+exp(-(y[3]+v7)/v8)));
// Imi current
double miinf2=1/(1+exp((y[3]-vhalf_mi)/k_mi));

double s_inf1 = 0.5*(1+tanh((y[0]-V_half)/V_slope));
double s_inf2 = 0.5*(1+tanh((y[3]-V_half)/V_slope));
//double Isyn1 = gsyn*y[7]*(Esyn-y[0]);	
//double Isyn2 = gsyn*y[6]*(Esyn-y[3]);	
double Isyn1 = gsyn*s_inf2*(Esyn-y[0]);	
double Isyn2 = gsyn*s_inf1*(Esyn-y[3]);	

// ML neuron #1
f[0]= (1/C)*(gL*(EL-y[0]) + gbarCa*minf*(ECa-y[0]) + gbarK*y[1]*(EK-y[0]) + gbarh*y[2]*(Eh-y[0])+ gmi*miinf*(Emi-y[0]) + Iapp + Isyn1);// voltage
f[1]= lambdaw*(winf-y[1]); //Ik
f[2] = (h_inf-y[2])/tauh; //Ih
// ML neuron #2
f[3]= (1/C)*(gL*(EL-y[3]) + gbarCa*minf2*(ECa-y[3]) + gbarK*y[4]*(EK-y[3]) + gbarh*y[5]*(Eh-y[3])+ gmi*miinf2*(Emi-y[3]) + Iapp + Isyn2);// voltage
f[4]= lambdaw2*(winf2-y[4]); //Ik
f[5] = (h_inf2-y[5])/tauh2; //Ih

f[6]=(s_inf1-y[6])/(tau_s); //synaptic activation 1
f[7]=(s_inf2-y[7])/(tau_s); //synaptic activation 2


}
void euler()
{
 int i=0;
 system();
  while (i<Neq){y[i]+=dt*f[i]; i++;}
} 

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
mxArray *tempy = mxCreateDoubleMatrix(1,Neq,mxREAL);
y=mxGetPr(tempy);
mxArray *tempf = mxCreateDoubleMatrix(1,Neq,mxREAL);
f=mxGetPr(tempf);
// inputs   
TT   = *(mxGetPr(prhs[0]));
gmi = *(mxGetPr(prhs[1]));
gbarK = *(mxGetPr(prhs[2]));
gbarCa = *(mxGetPr(prhs[3]));
gbarh = *(mxGetPr(prhs[4]));
gL = *(mxGetPr(prhs[5]));
gsyn = *(mxGetPr(prhs[6]));
V_half = *(mxGetPr(prhs[7]));
Iapp = *(mxGetPr(prhs[8]));

int kk;
tgl=0;
dt=0.02;
N=(int)(TT/dt);
// outputs    
double *Vm1;  
double *Vm2;  
double *m1;
double *m2;
plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[1] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
plhs[3] = mxCreateDoubleMatrix(1,N,mxREAL);

Vm1 = mxGetPr(plhs[0]);
m1 = mxGetPr(plhs[1]);
Vm2 = mxGetPr(plhs[2]);
m2 = mxGetPr(plhs[3]);

y[0]=-50; y[1]=0.; y[2]=0; // initial conditions neuron #1
y[3]=50; y[4]=0.; y[5]=0; // initial conditions neuron #2
y[6]=0; y[7]=0;
	
kk=0;
for(int i=0; i<N; i++)
{
tgl+=dt;
inp=0;

euler();
Vm1[kk]=y[0];
m1[kk]=y[1];
Vm2[kk]=y[3];
m2[kk]=y[4];

if(kk>N) break;
kk++;
   }
   return;
} // end mexFunction()
