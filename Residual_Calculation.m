function [ Res ] = Residual_Calculation(LHS,U,RHS )
%CALCULATE the residual 
%LHS generally is the stiffness matrix
% U is the vector of unknown\
%RHS is the right hand side vector
Res=LHS*U-RHS;

end

