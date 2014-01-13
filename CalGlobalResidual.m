clear all
% Global parameter input ----------------------------
global s;
s=10^(-4);
%Residual saturation of the liquid phase (Marchand's paper1, Table2)
global S_lr;
S_lr=0.4;
%Residual saturation of the gas phase ( Marchand's paper1, Table2)
global S_gr;
S_gr=0.0;
%The parameter to calculate the capillary pressure
global PC_0;
PC_0=20;%bar
%Viscosity of the liquid phase(Marchand's paper1, Table2)
global Mu_l;
Mu_l=3.161*10^(-18);
%Viscosity of the gas phase(Marchand's paper1, Table2)
global Mu_g;
Mu_g=2.845*10^(-20);
%Van Genuchten parameters
global m;
global n;
n=1.49;
m=1-(1/n);
%Henry's value
global Hen;
Hen=0.765;%mol/bar/m^3
%Molar mass of component 2
global M_2;
M_2=0.01;%kg/mol
%Parameter to calculate N_L
global N_L_std;
N_L_std=100;
%Parameter to calculate X_m
global X_m_crit;
X_m_crit=0.8;
%Parameter to calculate X_M
global X_M_crit;
X_M_crit=0.85;
%Parameter to calculate N_L
global X_L_crit;
X_L_crit=0.99995;
%Parameter to calculate X_m and X_M
global P_crit;
P_crit=100*PC_0;
%Intrinsic Permeability
global K;
K=5*10^(-20);
%The diffusion coefficient of the component (1) in Gas/Liquid phase
global DG;
global DL;
DG=9.467;DL=9.467;%m^2/centrury
%RT= R*T
global RT;
RT=1./39.7;%mol/bar
global poro;
poro=0.15;
global Stp;
global tot;
Stp=1000;
tot=10^(-5);
%dBeta used to control the derivation of the (F(X+dB*X)-F(X))/dB*X
dB=10^(-7);
% Parameter for initial value control
P_I=100;
X_I=0.18;
% numerical control
max_lin_sol_iter = 1000;
% parameter to control the secondary variable calculation
global max_iter_S;
max_iter_S=1000;
%Parameter for FEM Discretization
theta = 0.5;  % implicit
%theta = 1.;% explicit 
% theta = 0.5;  % Crank-Nicolson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preprocessing
nodes      = dlmread('2p_borgeat_2d_tri_coarser_node.txt');

elements   = dlmread('2p_borgeat_2d_tri_coarser_elements.txt');
elements=elements+1;
% now the spacial dimension and element type
[nn,sdim]  = size(nodes);
[ne,dummy] = size(elements);

% -------------------------------------
%|                                     |
%|                                     |
%|                                     |
%|                                     |
% -------------------------------------

index1 = find(nodes(:,1) == 0);% find the X-cord is equal to zero, right hand side boundary
bc_rhs_nodes=nodes(index1,:);%get the coordinate of the nodes located on the right hand side 
[nbcn,dummy] = size( bc_rhs_nodes );
%X and P are the primary variables
%Where X is the total molar fraction of the light component
%      P is the mean pressure
bc_nodes_value_P=100*ones(nbcn,1);
bc_nodes_value_X=(10^(-4))*ones(nbcn,1);

P_ini=ones(nn,1);
P_ini=P_I*P_ini;

X_ini=ones(nn,1);
X_ini=X_I*X_ini;
U_0=[P_ini',X_ini']';
% Processing ----------------------------------

% time control
% time_steps = [0:10:400, 500:100:10000, 11000:1000:86400]';
% time_steps = [0:10:400, 500:100:10000]';
% time_steps = [0:100:86400]';
time_steps = [0:100:8640]';
[steps,dummy] = size(time_steps);
steps = steps - 1; 


% initialize previous and current solution matrix
% Primary unknown --------P mean pressure 
P_pre  = sparse(nn,1 ); 
P_cur  = sparse(nn,1 ); 
P_pre  = P_ini; 
%% Primary unknown --------X molar fraction of the light component
X_pre  = sparse(nn,1 ); 
X_cur  = sparse(nn,1 ); 
X_pre  = X_ini; 

%coupled unknown: U 
U_pre=[P_pre',X_pre']';
U_cur=[P_cur',X_cur']';



% storage space of unknowns after each time step
P_record      = sparse(nn, steps+1);
P_record(:,1) = P_pre;

X_record      = zeros(nn, steps+1);
X_record(:,1) = X_pre;

% THE COUPLED UNKNOWN
U_RECORD=[P_record',X_record']';

%Initialize the parameter of M which would be 
%M=zeros(2*nn,2*nn);
%H=zeros(2*nn,2*nn);
% cleaning of LHS and RHS
LHS11=sparse(nn,nn);
LHS21    = sparse(nn,nn);
LHS12=sparse(nn,nn);
LHS22=sparse(nn,nn);
LHS=sparse(2*nn,2*nn);
RHS1_U= sparse(nn,1 );
RHS1_D= sparse(nn,1 );
RHS2_U=sparse(nn,1);
RHS2_D= sparse(nn,1 );
RHS=sparse(2*nn,1);
check=zeros(nn,nn);
S_ini=zeros(ne,1 );%Saturation of the gas phase

%START the time iteration
for ti = 1 : steps;
    %Parameter for time step
    dt = (time_steps(ti+1) - time_steps(ti))/0.1; 
    %Initialize the secondary variables    
    %Calculate the value of Secondary variables
    %Loop for each cell
    %Define a vector to check the Saturation for each cell 
    S_check=zeros(ne,1);
    for ie = 1 : ne       
        % get the coordinates of connecting nodes
        sctr  = elements(ie,:);
        %o=sctr(1,1);
        coord = nodes(sctr,:) ; 
        %%%%%%%%%Interpolation in a tri element%%%%%%%%%%%%%%%%
        [nrows,ncols] = size(coord);
        I = ones(nrows,1);
        A = 0.5 * det([I coord]);
        %calculate the P and X in each cell which is the interpolated value 
        P_L=(1/A)*[A/3,A/3,A/3]*P_pre(sctr);
        X_L=(1/A)*[A/3,A/3,A/3]*X_pre(sctr);
        %%%%%%%%%%%%%%%%%%%%%%%Calculate the secondary variables%%%%%%%%%%%
        [X_m,X_M,N_G,N_L,S,PC,Flag]= IterCalSec(P_L,X_L);
        %%%%%%%%%%%%%%%%%%%%%%END Calculate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S_check(ie)=S;
    
        [Se_l,Se_g]=EffectSat(S);        
        [Kr_l,Kr_g]=RelatPermeab(Se_l,Se_g);
        Lamda_l=CalculateLamda(Kr_l,Mu_l);
        Lamda_g=CalculateLamda(Kr_g,Mu_g);
        
        N_12=poro*(S*N_G+(1-S)*N_L);
        N_22=(-1)*poro*(S*N_G+(1-S)*N_L);        
        % local mass matrix
        %|       |        | 
        %|       |  M12   | 
        %|----------------|
        %|       |  M22   |
        %|       |        |
        M12    =  (1.0/dt)*N_12*shapeshape_tri(  coord );
        
        M22=  (1.0/dt)*N_22*shapeshape_tri(  coord );
        %check(sctr,sctr)=check(sctr,sctr)+(1.0/dt)*N_12*shapeshape_tri(  coord );
        
        % local mass matrix
        %|       |        | 
        %|Disp11 |        | 
        %|----------------|
        %|Disp21 |        |
        %|       |        |
        % local advection matrix
        B_1=(Lamda_l*N_L*X_m)+(Lamda_g*N_G*X_M);
        B_2=(Lamda_l*N_L*(1-X_m))+(Lamda_g*N_G*(1-X_M));
        Disp_11 = dshapedshape_tri(coord, [B_1, 0.0; 0.0,B_1]);%Dt=1
        Disp_21 = dshapedshape_tri(coord, [B_2, 0.0; 0.0, B_2]);
        %LHS11(sctr,sctr)=LHS11(sctr,sctr)+theta*(Disp_11);
        %In order to convert to a symmetric matrix time differential to the
        % first equation
        LHS11(sctr,sctr)=LHS11(sctr,sctr)+theta*(Disp_11);
        LHS21(sctr,sctr)=LHS21(sctr,sctr)+theta*(Disp_21);
        LHS12(sctr,sctr)=LHS12(sctr,sctr)+(M12);
        LHS22(sctr,sctr)=LHS22(sctr,sctr)+(M22);
       
        
        %LHS_U=LHS*U_pre;
        %RHS1(sctr)=RHS1(sctr)+M12*X_pre(sctr);%+theta*Disp_11*P_pre(sctr);
        %RHS2(sctr)=RHS2(sctr)+M22*X_pre(sctr);%+theta*Disp_21*P_pre(sctr)
        %Calculate the RHS
        %Calculate the terms contains Capillary pressure on the right hand
        RHS_PC_G_LC=Calc_RHS_PC(S,PC,N_G,X_M,Lamda_g,1,coord);%CALIBRATE LC(light component) (1-S*(2-S))*PC*&&[N_G*X_M*K*Kr_g*poro/mu]&&
        RHS_PC_L_LC=Calc_RHS_PC(S,PC,N_L,X_m,Lamda_l,2,coord);%(S*(2-S))*PC*&&[N_L*X_m*K*Kr_l*poro/mu]&&
        RHS_PC_G_HC=Calc_RHS_PC(S,PC,N_G,1-X_M,Lamda_g,1,coord);%HC(heavy component) (1-S*(2-S))*PC*&&[N_G*(1-X_M)*K*Kr_g*poro/mu]&&
        RHS_PC_L_HC=Calc_RHS_PC(S,PC,N_L,1-X_m,Lamda_l,2,coord);%(S*(2-S))*PC*&&[N_L*(1-X_m)*K*Kr_l*poro/mu]&&
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%Calculate the terms contains molar fraction on the right hand
        RHS_X_G_LC=Calc_RHS_X(S,X_M,N_G,DG,coord );
        RHS_X_G_HC=Calc_RHS_X(S,(1-X_M),N_G,DG,coord );
        RHS_X_L_LC=Calc_RHS_X(S,(X_m),N_L,DL,coord );
        RHS_X_L_HC=Calc_RHS_X(S,(1-X_m),N_L,DL,coord );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %RHS1-------first part of right hand [(1/DT)*[A]-(1-THETA)*[D]]*U_pre 
        %RHS1_U is upper part RHS1_D is the down part
        RHS1_U(sctr)=RHS1_U(sctr)+M12*X_pre(sctr)-(1-theta)*Disp_11*P_pre(sctr);
        RHS1_D(sctr)=RHS1_D(sctr)+M22*X_pre(sctr)-(1-theta)*Disp_21*P_pre(sctr);

        %RHS2------second part of right hand [(1-theta)*Ft+theta*[F]]
        %RHS2_U IS THE UPPERPART RHS2_D is the down part
        RHS2_U(sctr)=RHS2_U(sctr)+(RHS_PC_G_LC+RHS_PC_L_LC-RHS_X_G_LC-RHS_X_L_LC);
        RHS2_D(sctr)=RHS2_D(sctr)+(RHS_PC_G_HC+RHS_PC_L_HC-RHS_X_G_HC-RHS_X_L_HC);
        %RHS1(sctr)=RHS1(sctr)+M12*X_pre(sctr)+(1-theta)*Disp_11*P_pre(sctr)-;
    end

    %CONSTRUCT THE GLOBAL STIFF MATRIX AND RHS
    LHS_COMP={LHS11 LHS12;LHS21 LHS22};
    LHS=cell2mat(LHS_COMP);    
    RHS1=Vector_Combine(RHS1_U,RHS1_D);
    
    RHS2=Vector_Combine(RHS2_U,RHS2_D);
    RHS=RHS1+RHS2; 
    %Imposing the Dirichlet boundary condition
    for ib = 1 : nbcn
        idx = index1(ib);
        % b(i) -= A(i,ii)*u_bar
        RHS = RHS - LHS(:,idx) * bc_nodes_value_P(ib);
        % A(ii,ii) -> xii
        x_i_i = LHS(idx,idx);        
        % A(ii,j ) =  0
        LHS(idx,:) = 0.0;
        % A(i ,ii) = 0
        LHS(:,idx) = 0.0;     
        % b(ii) = xii * u_bar
        RHS(idx) = x_i_i * bc_nodes_value_P(ib);
        % A(ii,ii) <- xii
        LHS(idx,idx) = x_i_i;
    end
    for ib=1:nbcn
        idx=index1(ib);
        RHS=RHS-LHS(:,idx+nn)*bc_nodes_value_X(ib);
        x_in_in=LHS(idx+nn,idx+nn);
        LHS(idx+nn,:) = 0.0;
        % A(i ,ii) = 0
        LHS(:,idx+nn) = 0.0;
        RHS(idx+nn) =x_in_in * bc_nodes_value_X(ib);
        LHS(idx+nn,idx+nn) = x_in_in;
    end
    %translate the sparse matrix to full matrix just for checking the rank
    %and the det
    LHS_dense=full(LHS);
    LHS_dense11=full(LHS11);
    LHS_dense12=full(LHS12);
    LHS_dense21=full(LHS21);
    LHS_dense22=full(LHS22);
    %change the format of U from [p1 p2 p3....pn X1 X2 X3...Xn]T to [p1 X1 p2 X2...pn Xn]T
    U1_pre=Matrix_Swap(U_pre,1,nn);
    %Calculate residual  
    Res_1=Residual_Calculation(LHS,U_pre,RHS);
    str = ['time step',num2str(ti)];
    disp(str);
end