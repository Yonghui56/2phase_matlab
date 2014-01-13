function [h_R]=CalcH(lamda,M,N,x)
%Input Lamda=K*Kr/mu M is the mol.mass x is the molar fraction 
%Used to calc one of the terms in right hand side 
h_R=lamda*M*(N^2)*x;