function [ Coeff_X] = Calc_Coeff_X(N,D,S)
%Caluc the right hand side of the coeff of the x MOLAR FRACTION
global poro;%is the porosity
Coeff_X=N*D*S*poro;



end

