clear;
rng(1); %Set seed for random generator

%%Generate data with Gaussian Kernel
%Generate random (p-dimensional) z's which are sources of variation
N = 50; %Number of data points 
p = 2; %Dimension of variation sources
I = eye(N);
One = ones(N);
dd = 0.02;
Z = rand(N,p); %N of p-dimenional points z

title_text1 = 'Originally generated 2-d data as variation source';
scatter_label2d(Z,title_text1,dd) %Plot scatter plot of Z

%Generate n-D embeded data at higher-dimensional space
n = 10; %Dimension of observed data
sigma_kernel = 0.3;
[X,K] = generate_GK_nd(Z,n,sigma_kernel);

%Plot PCA score to see if there is any nonlinear pattern in the generated
%observed data
pct = 0.9; %the percentage of threhold eigen-values
pc1 = 1; %The index of the first component to be plotted
pc2 = 2;
title_text2 = 'Principle components of';
[PC_ind,eig_values] = scatter_PCA_label2d(X,pc1,pc2,pct,title_text2,dd);

%Standardize each coordinates of X; in other words, for each column, substract means and divide by unbiased standard deviation
std_X = std(X,0,1); 
%X_std = (X-repmat(mean(X,1),N,1))/diag(std_X);
X_std = X;

%%Inverse KPCA
[V,D,U] = svd(X_std,'econ');
sigma = 0.0; %Variation of noise
proj_mat = 1/N*One+V*V'-n*sigma^2*I; %Matrix for projection of K to column space of X
K_ini = X_std*X_std'; %Initial K matrix (non-negative definite)
K_ini = (K_ini+K_ini')/2;
count = 0; %Count for number of iteration
h_K = zeros(N); %Initial Matrix of inverse kernel
proj_cent = I-1/N*One; %Projection Matrix for centering
esp = 1e-6; %Tolerance for convergence
K_est = K_ini; 

sigma_alg = 0.3;
while 1
    %Projection
    K_tilde = proj_mat*K_est*proj_mat; %non-negative definite
    K_tilde = (K_tilde+K_tilde')/2;
    %Inverse function of kernel
    Z_est_new = IGaussian_Kernel(K_tilde,sigma_alg,p);
    %Condition of stopping loop because of convergence of estimation of Z
    count = count + 1;
    if count > 1
        diff_Z_est = norm(Z_est_new - Z_est);
        display(diff_Z_est);
        if diff_Z_est < esp
            break;
        end
    end
    Z_est = Z_est_new;
    %Update estimation of K based on new estimation of Z
    K_est = Gaussian_Kernel(Z_est,sigma_alg);
    if count == 1
        break;
    end
end

title_text3 = ['Final estimated variation source based on data in feature space'];
scatter_label2d(Z_est,title_text3,dd) %Plot scatter plot of estimation of Z