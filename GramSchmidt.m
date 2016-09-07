function [B] = GramSchmidt(A)
% Perform Gram-Schmidt orthogonalization onto column vectors of A
[m,n] = size(A);
[U, jb] = rref(A);
x = length(jb); 
C = zeros(m,x);
for i = 1:x
   C(:,i)= A(:,(jb(i))); 
end
B=C;
for i = 2:x
    for j = 1:i-1
        B(:,i) = C(:,i)- dot(C(:,i),B(:,j))/dot(B(:,j),B(:,j))* B(:,j) ;
    end
end
end