% Try drtoolbox functions

% Add the directory of drtoolbox to the search path of MATLAB and save it
addpath(genpath('/Applications/MATLAB_R2016a.app/toolbox/drtoolbox'))
savepath
%% An example for using the function
% Generate Helix data and draw a scatter plot
% Labels are just the colors used to label each data point: When C is a 
% vector the same number of rows as X (rows of X are instances), the values 
% in C are linearly mapped to the colors in the current colormap. When C is 
% a nrow(X)-by-3 matrix, the values in C specify the colors of the markers 
% as RGB values. C can also be a color string.
[X, labels] = generate_data('helix', 2000); 
figure, scatter3(X(:,1), X(:,2), X(:,3), 5, labels); title('Original dataset'), drawnow
% Estimate the intrinsic number of dimensions
% Use the MLE method
no_dims = round(intrinsic_dim(X, 'MLE'));
disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
% Use PCA method to do dimension reduction and plot 2-D PCA scores
[mappedX, mapping] = compute_mapping(X, 'PCA', no_dims);
figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels); title('Result of PCA');
% Use Laplacian Eigenmaps to do dimension reduction and plot 2-D mapped
% points. The number of nearest neighbor, k, is set to 7.
[mappedX, mapping] = compute_mapping(X, 'Laplacian', no_dims, 7);
figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels(mapping.conn_comp)); title('Result of Laplacian Eigenmaps'); drawnow
% Use
[mappedX, mapping] = compute_mapping(X, 'Isomap', no_dims);
figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels(mapping.conn_comp)); title('Result of Isomap'); drawnow
[mappedX, mapping] = compute_mapping(X, 'PCA', 2);

%% Try functions in drtoolbox on generated Gaussian profiles
options = ini_options();
%Generate random (p-dimensional) z's which are sources of variation
N = options.N; %Number of data points 
p = options.p; %Dimension of variation sources
dd = options.dd;
psize = options.psize;
rng(2); % Set seed for random generator
Z = rand(N,p); %N of p-dimenional points z
    
title_text1 = 'Originally generated 2-d data as variation source';
scat_func_tag = 1;
scatter_label2d_func = @scatter_label2d;
if scat_func_tag == 1
    tag = 0; %tag for showing index of data in 2d scatter plot (0: not showing)
    color = [0.7*ones(N,1),Z(:,1)/max(Z(:,1)),Z(:,2)/max(Z(:,2))];
    scatter_label2d_func(Z,title_text1,dd,tag,psize,color) %Plot scatter plot of Z
elseif scat_func_tag == 2
    tag = 1;
    color = [0.7*ones(N,1),Z(:,1)/max(Z(:,1)),0.3*ones(N,1)];
    a = int16(Z(:,2)*100); b = num2str(a); c = cellstr(b); %Labels
    scatter_label2d_func(Z,title_text1,dd,tag,psize,color,c) %Plot scatter plot of Z
end

% Generated data set
%Generate N Gaussian profile images (2-D) as high-dimensional observations
%sigma_data: the sigma of those Gaussian profiles
%l: the largest index of pixel in one side of images of those Gaussian
sigma_data = options.sigma_data;
l = options.l;
n = (l+1)^2;
X = generate_Gaussian_profile(Z*l,N,sigma_data,l);

%Plot PCA score to see if there is any nonlinear pattern in the generated
%observed data
pct = options.pct; %the percentage of threhold eigen-values
pc1 = options.pc1; %The index of the first component to be plotted
pc2 = options.pc2;
pc3 = options.pc3;
az = options.az;
el = options.el;
title_text2 = 'Principle components of';
color_3D = [0.7*ones(N,1),Z(:,1)/max(Z(:,1)),Z(:,2)/max(Z(:,2))];
%PCA on original observations
[PC_ind,eig_values] = scatter_PCA_3d(X,pc1,pc2,pc3,pct,title_text2,color_3D,psize,az,el);
%In svd (singular value decomposition, the singular values are sorted
%in non-increasing order.
[V,D,U] = svd((eye(N)-ones(N)/N)*X,'econ');
%Principle components of X (PCA scores)
PC_X = V*D(:,1:length(PC_ind));
%Plot 2-D PCA scores of X from my own code
title_text3 = ['Result of PCA by my own code'];
if scat_func_tag == 1
    scatter_label2d_func(PC_X(:,1:p),title_text3,dd,tag,psize,color) 
elseif scat_func_tag == 2
    scatter_label2d_func(PC_X(:,1:p),title_text3,dd,tag,psize,color,c) 
end
saveas(gca,[options.cwd,['PCA-my-own-code']],'jpg');
saveas(gca,[options.cwd,['PCA-my-own-code']],'fig');
close(gcf);

%Standardize each coordinates of X; in other words, for each column,
%substract means and divide by unbiased standard deviation of that column
std_X = std(X,0,1); 
X_std = (X-repmat(mean(X,1),N,1))/diag(std_X);

%PCA on standardized observations
title_text3 = 'Standardized version: Principle components of';
[PC_std_ind,eig_values_std] = scatter_PCA_3d(X_std,pc1,pc2,pc3,pct,title_text3,color_3D,psize,az,el);
%In svd (singular value decomposition, the singular values are sorted
%in non-increasing order.
[V,D,U] = svd(X_std,'econ');
%Principle components of X (PCA scores)
PC_X_std = V*D(:,1:length(PC_std_ind));
%Plot 2-D PCA scores of X_std from my own code
title_text3 = ['Result of PCA by my own code (standardized version)'];
if scat_func_tag == 1
    scatter_label2d_func(PC_X_std(:,1:p),title_text3,dd,tag,psize,color) 
elseif scat_func_tag == 2
    scatter_label2d_func(PC_X_std(:,1:p),title_text3,dd,tag,psize,color,c) 
end
saveas(gca,[options.cwd,['PCA(sv)-my-own-code']],'jpg');
saveas(gca,[options.cwd,['PCA(sv)-my-own-code']],'fig');
close(gcf);

% Running functions from drtoolbox
methods_name_list = {'PCA';'Isomap';'Laplacian';'HessianLLE';'LLE';'LDA';'MDS';'LandmarkIsomap';'LTSA';'KernelPCA';'tSNE';'LLTSA'};
no_dims = round(intrinsic_dim(X, 'MLE'));
for i = 1:length(methods_name_list)
    [mappedX, mapping] = compute_mapping(X, methods_name_list{i}, no_dims);
    figure, scatter(mappedX(:,1), mappedX(:,2), 5, color_3D); title(['Result of ', methods_name_list{i}]);
    saveas(gca,[options.cwd,methods_name_list{i}],'jpg');
    saveas(gca,[options.cwd,methods_name_list{i}],'fig');
    close(gcf);
end
for i = 1:length(methods_name_list)
    [mappedX, mapping] = compute_mapping(X_std, methods_name_list{i}, no_dims);
    figure, scatter(mappedX(:,1), mappedX(:,2), 5, color_3D); title(['Result of ', methods_name_list{i},' (standardized version)']);
    saveas(gca,[options.cwd,[methods_name_list{i},'(sv)']],'jpg');
    saveas(gca,[options.cwd,[methods_name_list{i},'(sv)']],'fig');
    close(gcf);
end
