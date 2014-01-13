% TODO: explain...
%global X_m_p_0; 
% TODO: explain...
%global X_M_p_1
% TODO: explain...

% residual value to control the iteration
global tot;
% parameter to calc the N_G
global s;
% parameter to calc the N_L (molar density)
global N_L_std;
% residual saturation of the liquid phase 
global S_lr;
% residual saturation of the gas phase 
global S_gr;
% parameter R*T
global RT;
% Parameter to calculate X_M.
global P_crit;
%Parameter to calc the capillary pressure
global PC_0;
%maximum iteration steps
global max_iter_S;

%parameter used for calculating the X_m X_M
global X_m_crit;
global X_M_crit;
%Parameter of Von model
global m;
global n;



n=1.49;
m=1-(1/n);
X_M_crit=0.85;
X_m_crit=0.8;
tot=10^(-5);
s=10^(-4);
N_L_std=85*10^3;
S_lr=0.0;
S_gr=0.0;
RT=1./39.7;
PC_0=20;
P_crit=100*PC_0;
max_iter_S=1000;
P = 7*1; 
X = 1*0.003;
[X_m,X_M,N_G,N_L,S,Flag] = IterCalSec(P,X); 