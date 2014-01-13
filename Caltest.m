function [S,G]=Caltest(Se_l,m,n,Se_g)
[Kr_l,Kr_g]=RelatPermeab(Se_l,Se_g,m);
G=Kr_l+Kr_g;
S=CapillaryP(Se_l,m,n);



