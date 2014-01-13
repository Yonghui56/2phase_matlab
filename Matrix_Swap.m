function [ B ] = Matrix_Swap(A,dimension,node )
B=A;
for i=1:node
    if dimension==1% which means A is the one dimension vector  
        
        B(2*i-1)=A(i);
        B(2*i)=A(i+node);    
    elseif dimension==2 %which means A is a two D matrix and we need to swap based on columns
        
        B(:,2*i-1)=A(:,i);
        B(:,2*i)=A(:,i+node);
    elseif dimension==3 % which means A is two D matrix and we need to swap based on rows 
        
        B(2*i-1,:)=A(i,:);
        B(2*i,:)=A(i+node,:);
    end
end





