function[Se_l,Se_g] = EffectSat(S)
% Function : calc the effective satutation of the gas phas and the liquid
% phase
% first evaluate the relation of the S and the residual S 
%Input : S saturation
%output Se_l is the effective saturation of the liquid phase
%output Se_g is the effective saturation of the gas phase

global S_lr;%residual saturation (liquid)
global S_gr;%residual saturation (gas)

if S>1-S_lr
    Se_l=0;
    Se_g=(S-S_lr)/(1-S_gr-S_lr);
elseif S<S_gr
        Se_l=(1-S-S_lr)/(1-S_gr-S_lr);
        Se_g=0;
else
    Se_l=(1-S-S_lr)/(1-S_gr-S_lr);
    Se_g=(S-S_gr)/(1-S_gr-S_lr);
end
    
        
