function [PC_ind,eig_values]=scatter_PCA_3d(X,pc1,pc2,pc3,pct,title_text,color_3D,psize,az,el)
%Plot 3d scatter plot of principle components of X
%Input:
%X: The columns are different coordinates of data; and rows are samples of data
%pc1: The index of the first component to be plotted
%pc2: The index of the second component to be plotted
%pc3: The index of the third component to be plotted
%pct: the percentage of threhold eigen-values
%title_text: The title text used in plot
%color_3D: color map for 3D plot (inheritate from 2D color map to show the
%   mapping from variation sources to observations
%psize: size of marker of data points.
%az: azmuthal angle of view
%el: elevation angle of view
%Output: A figure of scatter plot of PCA score

N = size(X,1);
%Singular Value Decomposition
[V,D,U] = svd((eye(N)-ones(N)/N)*X/sqrt(N-1),'econ');
%display(diag(D).^2);

%Find threshold index of singular value for which the sum of the corresponding
%eigenvalues, larger or equal to it, accounts for at least pct percentage
%of total variation.
th_ind = sing_th_ind(diag(D),pct);

PC_ind = (1:th_ind)
eig_values = diag(D).^2;

%Plot principle components with index pc1 and pc2 and pc3, according to the 
%printed eigenvalues
%(1) If there are some nonlinear curves in the plot in any pair of index, 
%then there will be some nonlinear manifold inside the data
%(2) If data are randomly scattering in all pairs of principle components,
%then data are pretty linearly in those principle components.
PCA_score = V*D; %PCA score of data
% display('The eigenvectors of covariance of X')
% U(:,1)
% norm(U(:,1))
% U(:,2)
% norm(U(:,2))
% U(:,3)
title_text = [title_text,' ',num2str(pc1),' , ',num2str(pc2),' and ',num2str(pc3),' PC'];
figure();
hold on;
subplot(1,2,1);
scatter_3d(PCA_score(:,[pc1,pc2,pc3]),title_text,color_3D,psize);
xlabel(num2str(pc1));
ylabel(num2str(pc2));
zlabel(num2str(pc3));
view(az,el);
subplot(1,2,2);

plot(1:2*th_ind,eig_values(1:2*th_ind),'*-');
axis square;
title(['PCA trend: threshold index = ',num2str(th_ind)]);
xlabel('Index of eigenvalues');
ylabel('Eigenvalues');
end