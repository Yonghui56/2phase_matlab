function [ PC ] = CalcCapillaryP( S )
%formulation to calculate the capillary pressure according to the Van
%Genuchten model
%S [0,1) when S is getting closed to 1 the value would be inf
%this function is used for the Function CapillaryP(S)
global m;
global n;
global PC_0;
[Se_l,~] = EffectSat(S);
%PC=(PC_0/10)*(Se_l^(-1/m)-1)^(1/n);
PC=PC_0*(Se_l^(-1/m)-1)^(1/n);
%PC=PC_0*S;

end

