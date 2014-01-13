%Calculate Saturation
%Input X primary variable Mean molar fraction of light component
%N_L is the molar density of the liquid phase
%N_G is the molar density of the gas phase
%X_M is the molar fraction of the light component in gas phase 
%X_m is the molar fraction of the light component in gas phase 
function[S]=CalcS(X,N_L,N_G,X_M,X_m)
S=(N_L)*(X-X_m)./(N_G*(X_M-X)+N_L*(X-X_m));