function [X_m,X_M,N_G,N_L,S,PC] = IterSatCal( P,X )
%This function is used to calculate the saturation as well as the other
%four secondary parameters based on the picard nonlinear Iteration
%   
global max_iter_S;
global tot;
%global S_gr;
%global S_lr;
%imposing the initial value of saturation
S_pi_Pre=0.015;
S_pi_Cur=S_pi_Pre+2*tot;
%str = ['step ', num2str(ip)];
%disp(str);
%record the iteration steps
i=1;
while abs(S_pi_Pre-S_pi_Cur)>=tot && i<max_iter_S
    S_pi_Pre=S_pi_Cur;
    
%     if S_pi_Pre<=1-S_lr && S_pi_Pre>=S_gr
%         S_bar=CalcS_bar( S_pi_Pre );
%         PC=CapillaryP(S_bar);
%     elseif S_pi_Pre<S_gr
%         S_gr_bar=CalcS_bar( S_gr );
%         PC=CapillaryP(S_gr_bar)+CalcDevPC(S_gr_bar)*(S_gr-S_pi_Pre);
%     elseif S_pi_Pre>1-S_lr
%         S_lr_bar=CalcS_bar( 1-S_lr );
%         PC=CapillaryP(S_lr_bar)+CalcDevPC(S_lr_bar)*(S_pi_Pre-1+S_lr);
%     end
    PC=CapillaryP(S_pi_Pre);
    [PG,PL]=CalcP(P,S_pi_Pre,PC);
        
    X_m= CalcX_m_Liq(PL,PC);
    %X_m=max(0,X_m_crit*(PG/(PG+P_crit)));
    X_M=CalculateX_M_gas(PG);
    N_G=CalcN_G(PG);
    N_L=CalcN_L(PL,X_m);
    S_pi_Cur=CalcS(X,N_L,N_G,X_M,X_m); 
    i=i+1;% Steps of the iteration
end
S= S_pi_Cur;   
[PG,PL]=CalcP(P,S_pi_Pre,PC);       
X_m= CalcX_m_Liq(PL,PC);
%X_m=max(0,X_m_crit*(PG/(PG+P_crit)));
X_M=CalculateX_M_gas(PG);
N_G=CalcN_G(PG);
N_L=CalcN_L(PL,X_m);        


