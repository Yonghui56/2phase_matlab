function [ S_bar ] = CalcS_bar( S )
%used to modify the capillary pressure 
% introduce a intermedia value(S_bar) to replace the saturation
%global parameter 
global S_gr;% saturation residual in the gas phase 
global S_lr;% saturation residual in the liquid phase 
eps=0.59;% is a parameter to control the value
S_bar=S_gr+(1-eps)*(S-S_gr)+0.5*eps*(1-S_gr-S_lr);

end

