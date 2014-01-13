%function calculate the molar fraction of the light component in gas phase
function[X_M]=CalculateX_M_gas(PG)
%Input PG is the pressure in gas phase
%Output X_M is the molar fraction of the light component in gas phase
% formulation comes from Marchand's paper
global X_M_crit
global P_crit;
[X_M]=min(X_M_crit,X_M_crit*(PG/P_crit)^(1/2));