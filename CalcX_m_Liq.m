%This function is used to calculate the molar fraction of light componentthe in liquid phase 

function[X_m]=CalcX_m_Liq(PL,PC)
%Input PL is the pressure in Liquid phase
%Input PC is the capillary pressure
%output X_m is the molar fraction of the liquid phase (light component)
global X_m_crit;
global P_crit;
X_m=max(0,X_m_crit*((PL+PC)/(PL+PC+P_crit)));