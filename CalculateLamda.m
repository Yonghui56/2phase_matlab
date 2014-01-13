%USE TO CALC K*Kr/mu
function [Lamda]=CalculateLamda(Kr,mu)
%absolute permeability
%Kr is relative permeability
%mu is viscosity
global K;
Lamda=K*Kr/mu;