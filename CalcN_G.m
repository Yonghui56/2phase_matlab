%Funtion calculate the molar density of the gas phase
function[N_G]=CalcN_G(PG)
%Input PG is the pressure in gas phase
%output N_G the molar density of the gas phase
global RT;
N_G=PG/RT;