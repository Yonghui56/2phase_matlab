function [RHS_PC] = Calc_RHS_PC(S,PC,N,X,Lamda,flag,coord )
% Function used to calc the Coeff of Capillary pressure on the right hand
% side
%Input S saturation PC capillary pressure,N molar density lamda K*Kr/mu
%flag control the formulation of PC
if 1==flag
    y_PC=(1-S*(2-S))*PC;
elseif 2==flag
    y_PC=S*(2-S)*PC;
end
Coeff_PC=Calc_Coeff_PC(N,X,Lamda);
RHS_PC=dshapedshape_tri(coord, [Coeff_PC, 0.0; 0.0, Coeff_PC])*[y_PC;y_PC;y_PC];


end

