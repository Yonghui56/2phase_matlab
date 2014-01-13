function [ devPC ] = CalcDevPC(S_eps_bar,S )
%Derivative of the function of CalcCapillaryP
%S is saturation
global eps;

devPC=(CalcCapillaryP(S_eps_bar)-CalcCapillaryP(S))/eps;


end

