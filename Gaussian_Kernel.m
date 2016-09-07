function K = Gaussian_Kernel(Z,sigma_kernel)
% Gaussian Kernel
% Input:
% Z: matrix of input data points with rows as points
% sigma_kernel: parameter for Gaussian Kernel
% Output:
% K: output Gaussian Kernel based on input data points

N = size(Z,1); %Number of data points
%Apply Gaussian Kernel to distance of two points from Z
K = zeros(N);
for i = 1:N
    for j = i:N
        diff_ij = Z(i,:)-Z(j,:);
        K(i,j) = exp(-diff_ij*diff_ij'/2/sigma_kernel^2);
        if i~=j
            K(j,i) = K(i,j);
        end
    end
end