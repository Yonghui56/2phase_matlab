function [ RHS_X] = Calc_RHS_X(S,X,N,D,coord )
%UNTITLED Summary of this function goes here
% input S saturation X molar fraction N molar density D diffusion coeff
Coeff_X=Calc_Coeff_X(N,D,S);%Calc the coeff of X
RHS_X=dshapedshape_tri(coord, [Coeff_X, 0.0; 0.0, Coeff_X])*[X;X;X];


end

