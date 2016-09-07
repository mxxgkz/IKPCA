clear;
rng(2); %Set seed for random generator
cwd = '/Users/kungangzhang/Documents/OneDrive/Northwestern/Study/Courses/Independent Study/20160713-Test-alg-on-new-generated-data/figures/';

%%Method1: generate the linear embedded data and nonlinear embedded data 
%%(the assumed model) apply PCA to them to see the principle components.

%%Generate data with Gaussian Kernel
%Generate random (p-dimensional) z's which are sources of variation
N = 1000; %Number of data points 
p = 2; %Dimension of variation sources
I = eye(N);
One = ones(N);
dd = 0.02;
Z = rand(N,p); %N of p-dimenional points z

tag = 0;

title_text1 = 'Originally generated 2-d data as variation sources';
scatter_label2d(Z,title_text1,dd,tag) %Plot scatter plot of Z

% %Generate n-D linearly embedded data at higher-dimensional space to test
% %the PCA and plot function
% n = 10; %Dimension of observed data
% X_base = rand(p,n);
% X_base_orth = (GramSchmidt(X_base'))';
% X_base_orth(1,:) = X_base_orth(1,:)/norm(X_base_orth(1,:));
% X_base_orth(2,:) = X_base_orth(2,:)/norm(X_base_orth(2,:));
% X_l = Z*X_base_orth;


%Generate n-Dim embedded data at higher-dimensional space
% n = 10; %Dimension of observed data
% sigma_kernel = 0.3; 
% [X,K] = generate_GK_nd(Z,n,sigma_kernel);

%Generate 2-D images of Gaussian profiles
%sigma_data: the sigma of those Gaussian profiles
%l: the largest index of pixel in one side of images of those Gaussian
sigma_data = 10;
l = 30;
X = generate_Gaussian_profile(Z*l,N,sigma_data,l);

%Plot PCA score to see if there is any nonlinear pattern in the generated
%observed data
pct = 0.9; %the percentage of threhold eigen-values
pc1 = 1; %The index of the first component to be plotted
pc2 = 2;
pc3 = 3;
title_text2 = 'Principle components of';
az = 30;
el = 25;
[PC_ind,eig_values] = scatter_PCA_3d(X,pc1,pc2,pc3,pct,title_text2,az,el);

% view(az,el);
% saveas(gcf,[cwd,'0'],'jpg');
% saveas(gcf,[cwd,'0'],'fig');

% title_text_l = 'Linear-embedding-Principle components of';
% [PC_ind_l,eig_values_l] = scatter_PCA_3d(X_l,pc1,pc2,pc3,pct,title_text_l);

%Standardize each coordinates of X; in other words, for each column, substract means and divide by unbiased standard deviation
std_X = std(X,0,1); 
%X_std = (X-repmat(mean(X,1),N,1))/diag(std_X);
X_std = X;