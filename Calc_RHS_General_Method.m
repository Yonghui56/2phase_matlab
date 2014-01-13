function [bt ] = Calc_RHS_General_Method(f, coord )
%
A = 0.5 * det([I coord]);% the area of the triangle cell
mid1 = (coord(1,:)+coord(2,:))/2;%Calc the coordinate of the midpoint of each edge
mid2 =(coord(1,:)+coord(3,:))/2;
mid3 = (coord(2,:)+coord(3,:))/2;
bt1 = A.*(f(mid2)+f(mid3))/6;
bt2 = A.*(f(mid3)+f(mid1))/6;
bt3 = A.*(f(mid1)+f(mid2))/6;
bt=[bt1;bt2;bt3];

end

