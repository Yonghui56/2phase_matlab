%CALCULATE the B_L 
%Input PL pressure in the liquid phase
%Input X_m molar fraction of the liquid phase 
%s is the parameter for example:s=10^(-4)
%N_L_std  is the parameter to calculate N_L
function[N_L]=CalcN_L(PL,X_m)
global s;
global N_L_std
B_L=min(1,1/(1+s*PL));
N_L=N_L_std/B_L/(1-X_m);