function [X,K] = generate_GK_nd(Z,n,sigma_kernel)
% Generate n-D embeded data at higher-dimensional space
% Input:
% Z: Data points of variation sources (p-dimensional)
% n: Dimension of higher-dimensional space where data Z are embeded
% Output:
% X: Data points of Z embeded in higher-dimensional space
% K: output Gaussian Kernel based on input data points

rng(1); %Step seed for random generator
N = size(Z,1); %Number of data points
%Generate random matrix
A = rand(n,N);
while (rank(A)<10) %If A doesn't have full rank, regenerate A. The reason is we assume X has full column rank.
    A = rand(n,N);
end

A_tilde = (GramSchmidt(A'))'; %Perform Gram-Schmidt orthogonalization on to rows of A???

%Use Gaussian Kernel to generate observed data X
K = Gaussian_Kernel(Z,sigma_kernel);
X = (A_tilde*K)'; %Linear combination of columns of K