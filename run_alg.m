function run_alg(rngn, scatter_label2d_func, scat_func_tag, options)
%%
% This function is to run the algorthm on Notes-20160520 on a data set of
% Gaussian profile. The goal is to look the mapping between variation
% sources and estimated variation sources, so that the performance of this
% algorithm can be evaluated. The output of this function is mainly
% scatter plots showing images. Two possible plot functions can be used:
% one is using 2-dimensional color map to code the scatter plot; the other
% use color and text to code the scatter plot.
%
% Input:
%   rngn: the random number generator
%   scatter_label2d_func: the name of function of scatter plot. two
%       possible choice: scatter_label2d (use 2-D color map);
%       scatter_label2d_1 (use color and text)
%   scat_func_tag: tag for two possible functions of scatter plot. 
%       scatter_label2d (1);
%       scatter_label2d_1 (2)
%   options: include all parameters in the algorithm
%%
    rng(rngn); %Set seed for random generator

    %Generate random (p-dimensional) z's which are sources of variation
    N = options.N; %Number of data points 
    p = options.p; %Dimension of variation sources
    dd = options.dd;
    psize = options.psize;
    
    I = eye(N);
    One = ones(N);
    Z = rand(N,p); %N of p-dimenional points z
    
    title_text1 = 'Originally generated 2-d data as variation source';
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
    % Generated data set 1
    % %Generate n-D embeded data at higher-dimensional space
    % n = 10; %Dimension of observed data
    % sigma_kernel = 0.3;
    % [X,K] = generate_GK_nd(Z,n,sigma_kernel);

    % Generated data set 2
    %Generate N Gaussian profile images (2-D) as high-dimensional observations
    %sigma_data: the sigma of those Gaussian profiles
    %l: the largest index of pixel in one side of images of those Gaussian
    sigma_data = options.sigma_data;
    l = options.l;
    n = (l+1)^2;
    X = generate_Gaussian_profile(Z*l,N,sigma_data,l);

    %Old plot
    % %Plot PCA score to see if there is any nonlinear pattern in the generated
    % %observed data
    % pct = 0.9; %the percentage of threhold eigen-values
    % pc1 = 1; %The index of the first component to be plotted
    % pc2 = 2;
    % title_text2 = ['Principle components of ',num2str(pc1),' and ',num2str(pc2),' PC'];
    % PC_ind = scatter_PCA_label2d(X,pc1,pc2,pct,title_text2,dd);

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
    
    %Standardize each coordinates of X; in other words, for each column,
    %substract means and divide by unbiased standard deviation of that column
    std_X = std(X,0,1); 
    X_std = (X-repmat(mean(X,1),N,1))/diag(std_X);
    %X_std = X;
    
    %PCA on standardized observations
    title_text3 = 'Standardized version: Principle components of';
    [PC_std_ind,eig_values_std] = scatter_PCA_3d(X_std,pc1,pc2,pc3,pct,title_text3,color_3D,psize,az,el);

    %%Inverse KPCA
    %In svd (singular value decomposition, the singular values are sorted
    %in non-increasing order.
    [V,D,U] = svd(X_std/sqrt(N-1),'econ');
    %Principle components of X (PCA scores)
    PC_X = V*D(:,1:length(PC_std_ind));
    %The standardized observations doesn't have full column rank, so we
    %should pick out those columns of V that correspond to non-zero
    %singular values.
    rank_X_std = length(PC_std_ind);
    V = V(:,1:rank_X_std);
    

    
    sigma_nois = options.sigma_nois; %Standard Variation of noise
    proj_mat = 1/N*One+V*V'-n*sigma_nois^2*I; %Matrix for projection of K to column space of X
    K_ini = X_std*X_std'; %Initial K matrix (non-negative definite)
    K_ini = (K_ini+K_ini')/2;
    count = 0; %Count for number of iteration
    h_K = zeros(N); %Initial Matrix of inverse kernel
    proj_cent = I-1/N*One; %Projection Matrix for centering
    esp = options.esp; %Tolerance for convergence
    K_est = K_ini; 

    
    sigma_alg = options.sigma_alg;
    while 1
        %Projection
        K_tilde = proj_mat*K_est*proj_mat; %non-negative definite
        K_tilde = (K_tilde+K_tilde')/2;
        %Inverse function of kernel
        [Z_est_new,lambda_est_new] = IGaussian_Kernel(K_tilde,sigma_alg,p);
        [Z_est_plot,lambda_est_plot] = IGaussian_Kernel(K_tilde,sigma_alg,4);
        contour_Z(Z(:,1),Z(:,2),real(Z_est_plot),lambda_est_plot,color_3D,psize);
        %Condition of stopping loop because of convergence of estimation of Z
        count = count + 1;
        if count > 1
            diff_Z_est = norm(Z_est_new(:,1:p) - Z_est(:,1:p));
            display(diff_Z_est);
            if diff_Z_est < esp
                break;
            end
        end
        saveas(gca,[options.cwd,['8-',num2str(count-1),'1']],'jpg');
        saveas(gca,[options.cwd,['8-',num2str(count-1),'1']],'fig');
        close(gcf);
        Z_est = Z_est_new;
        %Plug in the initial Z to see effect
        if count == 1
%             Z_est = PC_X;
            Z_est = Z;
        end
        
        %Update estimation of K based on new estimation of Z
        K_est = Gaussian_Kernel(Z_est,sigma_alg);
        
        %Plot the scattering plot of estimated Z
        title_text3 = [num2str(count-1),' estimated variation source based on data in feature space'];
        if scat_func_tag == 1
            scatter_label2d_func(Z_est(:,1:p),title_text3,dd,tag,psize,color) %Plot scatter plot of estimation of Z
        elseif scat_func_tag == 2
            scatter_label2d_func(Z_est(:,1:p),title_text3,dd,tag,psize,color,c) %Plot scatter plot of estimation of Z
        end
        saveas(gca,[options.cwd,['8-',num2str(count-1),'2']],'jpg');
        saveas(gca,[options.cwd,['8-',num2str(count-1),'2']],'fig');
        close(gcf);
        if count == options.max_iter
            break;
        end
        
        %pause;
        %K_est
    end

    title_text3 = ['Final estimated variation source based on data in feature space'];
    if scat_func_tag == 1
        scatter_label2d_func(Z_est,title_text3,dd,tag,psize,color) %Plot scatter plot of estimation of Z
    elseif scat_func_tag == 2
        scatter_label2d_func(Z_est,title_text3,dd,tag,psize,color,c) %Plot scatter plot of estimation of Z
    end
    
    count
end