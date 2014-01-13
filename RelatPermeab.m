% calculate the relative permeability 
% use the formulation of marchand's paper (2)(P696)
function[Kr_l,Kr_g] = RelatPermeab(Se_l,Se_g)
%Input Se_l effective saturation (liquid phase)
%Input Se_g effective saturation (gas phase)
%output Kr_l relative permeability (liquid phase)
%output Kr_g relative permeability (gas phase)
global m;


Kr_l=((Se_l)^(1/2))*((1-(1-Se_l^(1/m))^(m))^(2));
Kr_g=((Se_g)^(0.59))*((1-(1-Se_g^(1/m)))^(2*m));


