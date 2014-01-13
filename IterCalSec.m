%Use the value of the Primary variables to calculate the value of the secondary unknowns
%Five nonlinear equations and five unknowns 
%use the formulation of Marchand's paper
function[X_m,X_M,N_G,N_L,S,PC,Flag] = IterCalSec(P,X)

%%%%%%%%%%%%%%%%%  Input  %%%%%%%%%%%%%%%%%%%%%%%% 
%X is the primary unknowns 
%P is the primary unknowns
%X_m_p_0
%X_M_p_1
%Steps is the iteration steps
%%%%%%%%%%%%%%%%%  END Input  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Output  %%%%%%%%%%%%%%%%%%%%%%%%
%X_m is the component (1) in the liquid phase 
%X_M is the component (1) in the gas phase
%N_G is the Molar density of the gas phase
%N_L is the Molar density of the liquid phase
%S is the saturation of the gas phase 
% Flag to indicate the phase state 
%0---single liquid phase 1---single gas phase 2---two phase state 
%%%%%%%%%%%%%%%%%  END Output  %%%%%%%%%%%%%%%%%%%%
% residual value to control the iteration

% parameter to calc the N_G

% parameter to calc the N_L (molar density)
%global N_L_std;
% residual saturation of the liquid phase 
%global S_lr;
% residual saturation of the gas phase 
%global S_gr;
% parameter R*T
%global RT;
% PARAMETER FOR CALCULATE THE 
%global P_crit;
%parameter to control the iteration

%calculate the 
%X_m(p,0)---molar fraction of light component in Liquid phase when S =0  
%X_M(p,1)---molar fraction of light component in Gas phase when S =1

PC_s0=CapillaryP(0.0);%PC(S=0)

[~,PL_0]=CalcP(P,0,PC_s0);
PC_s1=CapillaryP(1.0);%PC(S=1)
[PG_1,~]=CalcP(P,1,PC_s1);
X_m_p_0=CalcX_m_Liq(PL_0,PC_s0);%(Liquid phase )
X_M_p_1=CalculateX_M_gas(PG_1);%(gas phase)
%Parameter to calc the capillary pressure
%global PC_0;
if X>X_m_p_0 && X<X_M_p_1
%PICARD ITERATION METHOD TO CALCULATE THE SECONDARY VARIABLES
    Flag=2;
    %calculate secondary variable 
    [X_m,X_M,N_G,N_L,S,PC ] = IterSatCal( P,X );
            
elseif X<=X_m_p_0 %only liquid phase
    S=0;
    Flag=0;    
    PL=P;
    N_G=0;
    X_m=X;
    X_M=0;    
    N_L=CalcN_L(PL,X_m);
    PC= PC_s0;       
elseif X>=X_M_p_1 % Only gas phase
    S=1;
    Flag=1;
    PG=P;
    X_M=X;
    X_m=0;
    N_G=CalcN_G(PG);
    N_L=0;
    PC = PC_s1;
end