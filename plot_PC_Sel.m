K=101;
S=linspace(0,1,K);
global S_lr;
global S_gr;
global m;
global n;
global PC_0;
PC_0=20;
n=1.49;
m=1-1/n;
S_lr=0.1;
S_gr=0.05;
global eps;
PC=zeros(K,1);
PC2=zeros(K,1);
S_lrJ=zeros(K,2);
eps=10e-8;
for i=1:K
    S_pi_Pre=S(i);
%     if S(i)>1
%         Y=1;
%     elseif S(i)<0
%         Y=0;
%     else
%         Y=S(i)*(2-S(i));
%     end
    if S_pi_Pre<=1-S_lr && S_pi_Pre>=S_gr
        S_bar=CalcS_bar( S_pi_Pre );
        PC(i)=CalcCapillaryP(S_bar);
    elseif S_pi_Pre==0
        PC(i)=0;
    elseif S_pi_Pre<S_gr && S_pi_Pre>0
        S_gr_bar=CalcS_bar( S_gr );
        S_gr_eps_bar=CalcS_bar(S_gr+eps);
        %PC(i)=CalcCapillaryP(S_gr_bar)+CalcDevPC(S_gr_eps_bar,S_gr_bar)*(S_gr-S_pi_Pre);
        PC(i)=(CalcCapillaryP(S_gr_bar)/S_gr)*S_pi_Pre;
    elseif S_pi_Pre>1-S_lr && S_pi_Pre<=1
        S_lr_bar=CalcS_bar( 1-S_lr );
        S_lr_eps_bar=CalcS_bar(1-S_lr+eps) ;
        PC(i)=CalcCapillaryP(S_lr_bar)+CalcDevPC(S_lr_eps_bar, S_lr_bar )*(S_pi_Pre-1+S_lr);
        S_lrJ(i,1)=S_lr_bar;
%     elseif S_pi_Pre>1
%         S_lr_bar=CalcS_bar( 1-S_lr );
%         S_lr_eps_bar=CalcS_bar( 1-S_lr+eps );
%         PC(i)=CapillaryP(S_lr_bar)+CalcDevPC(S_lr_eps_bar,S_lr_bar)*(S_lr);
%     elseif S_pi_Pre<0
%         PC(i)=0;
        
    end
end
for i=1:K
    PC2(i)=testCalc_PC(S(i));
end
plot(S,PC);
figure;
%plot(S,PC2);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 