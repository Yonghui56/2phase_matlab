function [PG,PL]=CalcP(P,S,PC)
%P is the mean pressure
%S is the saturation 
%PC is the capillary pressure 
%PG is the pressure in the gas phase 
%PL is the pressure in the liquid phase

if  S<=0%Abs(S)<10^(-6)S=0
    PG=0;
    PL=P;
elseif S>=1%Abs(S)> S=1
    PG=P;
    PL=0;
else
    % the formulation comes from the Marchand's paper(1) (P434)
    PG=P+(1-S*(2-S))*PC;
    PL=P-S*(2-S)*PC;
end
