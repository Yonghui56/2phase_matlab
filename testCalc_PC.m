	 
%Calculate the capillary pressure
%use the formulation of Marchand's paper(2)(P696)
function[PC] = testCalc_PC(S)
global PC_0;
global m;
global n;
eps=10^(-6);
if abs(S-0)<=eps% means S (gas phase) is equal to zero
    PC=0;
elseif abs(S-1)<=eps% means S (gas phase) is equal to one
    PC=PC_0; 
else
    [Se_l,~] = EffectSat(S);
    %PC=(PC_0/10)*(Se_l^(-1/m)-1)^(1/n);
    PC=PC_0*(Se_l^(-1/m)-1)^(1/n);
    %PC=PC_0*S;
end

