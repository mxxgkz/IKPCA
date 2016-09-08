function [options] = ini_options()
%%
% Initialize parameters in the run_alg.m
% options.N: sample size
% options.p: dimension of variation sources
% options.dd: the translation of text besides a data point in scatter plot
% options.psize: the size of marker in scatter plot
% options.sigma_alg: the sigma in algorithm (Gaussian Kernel)
% options.sigma_nois: the sigma of nois in algorithm when projecting kernel
%   matrix
% options.pct: the threshold percentage of variation in PCA on observations
%   X
% options.pc1,pc2,pc3: the index of three principle components
% options.az: azthmuthal angle of the view
% options.el: elevation angle of the view
% options.sigma_data: the sigma used to generate data set (images of
%   Gaussian profiles)
% options.l: length of images
% options.esp: the tolerance of norm of difference in Z between two 
%   consecutive iteration
% options.max_iter: maximum iteration number
%%
% Parameters for generating data
    options.N = 200; %Number of data points 
    options.l = 8;
    options.p = 2; %Dimension of variation sources
    options.sigma_data = 10;
    options.sigma_nois = 0.0; %Standard Variation of noise
    
% Parameters for running the algorithm    
    options.sigma_alg = 0.25;
    options.max_iter = 50;
    options.esp = 1e-4; %Tolerance for convergence

% Parameters for storing and plotting data
    options.cwd = '/Users/kungangzhang/Documents/OneDrive/Northwestern/Study/Courses/Independent Study/20160907-Implement-5-step-alg/figures/';
    options.pct = 0.9999; %the percentage of threhold eigen-values
    options.pc1 = 1; %The index of the first component to be plotted
    options.pc2 = 2;
    options.pc3 = 3;
    options.az = 30;
    options.el = 25;
    options.dd = 0.02;
    options.psize = 30;
        
end