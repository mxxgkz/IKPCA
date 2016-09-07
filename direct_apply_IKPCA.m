clear;
rng(1); %Set seed for random generator

%%Generate data with Gaussian Kernel
%Generate random (p-dimensional) z's which are sources of variation
N = 50; %Number of data points 
p = 2; %Dimension of variation sources
Z = rand(N,p); %N of p-dimenional points z
I = eye(N);
One = ones(N);
dd = 0.01;
Z = rand(N,p); %N of p-dimenional points z

scatter_label2d(Z,dd) %Plot scatter plot of Z

%Generate n-D embeded data at higher-dimensional space
n = 10; %Dimension of observed data
sigma_kernel = 0.1;
[X,K] = generate_GK_nd(Z,n,sigma_kernel);

%Apply Inverse KPCA directly to the original K, and the data are almost perfectly recovered.

[V,D,U] = svd(X,'econ');
sigma = 0.1; %Variation of noise

%%Update Z as an eigenvector problem
Z_est = IGaussian_Kernel(K,sigma_kernel,p);

scatter_label2d(Z_est,dd) %Plot scatter plot of estimation of Z