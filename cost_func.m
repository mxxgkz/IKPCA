function [f,g] = cost_func(Z,A,B,sigma)
% Calculate function value
[N,n] =size(Z);
f = 0;
K = Gaussian_Kernel(Z,sigma);
for i = 1:n
    f=f+B(i,:)-A(i,:)*
end

end