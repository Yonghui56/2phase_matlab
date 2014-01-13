%Calculate the capillary pressure 
%use the formulation of Marchand's paper(2)(P696)
function[PC] = CapillaryP(S)
%Input Se_l is the effective saturation of the liquid phase 
%m,n is the Von Genuchten parameter
%PC_0 is a parameter to calculate PC
%eps=10^(-6);
% if abs(S-0)<=eps% means S (gas phase) is equal to zero
%     PC=0;
% elseif abs(S-1)<=eps% means S (gas phase) is equal to one
%     PC=PC_0; 
% else
%     [Se_l,~] = EffectSat(S);
%     %PC=(PC_0/10)*(Se_l^(-1/m)-1)^(1/n);
%     PC=PC_0*(Se_l^(-1/m)-1)^(1/n);
%     %PC=PC_0*S;
% end
% extend the definition for S S could be >1 or <0
global S_gr;
global S_lr;
global eps;
% Implement on 8th Jan
%the low point used the slope of the point (Sgr)
% if S<=1-S_lr && S>=S_gr
%     S_bar=CalcS_bar( S );
%     PC=CalcCapillaryP(S_bar);
% elseif S<S_gr
%     S_gr_bar=CalcS_bar( S_gr );
% %   S_gr_eps_bar=CalcS_bar( S_gr+eps );
%     PC=CalcCapillaryP(S_gr_bar);%+CalcDevPC(S_gr_eps_bar,S_gr_bar)*(S_gr-S);
% elseif S>1-S_lr
%     S_lr_bar=CalcS_bar( 1-S_lr );        
%     PC=CalcCapillaryP(S_lr_bar);%+CalcDevPC(S_lr_eps_bar,S_lr_bar)*(S-1+S_lr));
%         
% end
%Implement on 9th Jan 
%The low point connect  point with (0,0) by a straight line 
    if S<=1-S_lr && S>=S_gr
        S_bar=CalcS_bar( S );
        PC=CalcCapillaryP(S_bar);
        %if S_gr is equal to zero
%         if S_gr==0
%             PC=0;
%         end
    elseif S==0
        PC=0;
    elseif S<S_gr && S>0
        S_gr_bar=CalcS_bar( S_gr );
        %S_gr_eps_bar=CalcS_bar(S_gr+eps);
        %PC(i)=CalcCapillaryP(S_gr_bar)+CalcDevPC(S_gr_eps_bar,S_gr_bar)*(S_gr-S);
        PC=(CalcCapillaryP(S_gr_bar)/S_gr)*S;
    elseif S>1-S_lr 
        S_lr_bar=CalcS_bar( 1-S_lr );
        S_lr_eps_bar=CalcS_bar(1-S_lr+eps) ;
        PC=CalcCapillaryP(S_lr_bar)+CalcDevPC(S_lr_eps_bar, S_lr_bar )*(S-1+S_lr);
        
    end