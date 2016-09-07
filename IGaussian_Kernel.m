function [Z,lambda] = IGaussian_Kernel(K,sigma_kernel,p)
% Inverse Gaussian Kernel
% Input:
% K: output Gaussian Kernel based on input data points
% sigma_kernel: parameter for Gaussian Kernel
% p: dimension of data points of Z
% Output:
% Z: matrix of input data points with rows as points

N = size(K,1); %Number of data points
%Inverse function of kernel, make sure h_K is a Hermitian
h_K = zeros(N);
for i = 1:N
    for j = i:N
        h_K(i,j) = -2*sigma_kernel^2*log(K(i,j));
        if i~=j
            h_K(j,i) = conj(h_K(i,j));
        end
    end
end
%MDS: recover estimation of Z (p-dimensional points)
I = eye(N);
One = ones(N);
proj_cent = I-1/N*One;
%diag(h_K)
h_K_new = -0.5*proj_cent*h_K*proj_cent;
h_K_new = (h_K_new + h_K_new')/2;
%diag(h_K_new)
%In eig funciton, the eigen-values are sorted in non-decreasing order
[eig_vectors,eig_values] = eig(h_K_new);
%display(diag(eig_values));
Z = (I(N:-1:N-(p-1),:)*sqrt(eig_values)*eig_vectors')';
lambda = diag(eig_values);
th_ind = sing_th_ind(sqrt(lambda(N:-1:1)),0.99);
lambda = lambda(N:-1:N-(p-1));
end