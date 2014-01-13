function [b_R]=CalcB(N,D,S)
%Used to calc one of the terms of the right hand side 
%N is the molar density D is the diffusion coeff S is the saturation
b_R=N*D*S;
