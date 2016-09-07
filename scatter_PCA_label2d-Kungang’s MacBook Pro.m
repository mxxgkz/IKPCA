function [PC_ind,eig_values] = scatter_PCA_label2d(X,pc1,pc2,pct,title_text,dd)
%Plot 2d scatter plot of principle components of X
%Input:
%X: The columns are different coordinates of data; and rows are samples of data
%pc1: The index of the first component to be plotted
%pc2: The index of the second component to be plotted
%pct: the percentage of threhold eigen-values
%title_text: The title text used in plot
%dd: the factor of shifted translation for labels
%Output: A figure of scatter plot of PCA score

N = size(X,1);
%Singular Value Decomposition
[V,D,U] = svd(X/sqrt(N-1),'econ');
display(diag(D));
sum_eig = sum(diag(D).^2); %Sum of all eigen-values
th_ind = 1; %The threshold index of principle components
tem_sum_eig = D(th_ind,th_ind)^2;

while th_ind <= size(D,1) && tem_sum_eig/sum_eig < pct
    th_ind = th_ind + 1;
    tem_sum_eig = tem_sum_eig + D(th_ind,th_ind)^2;
end
PC_ind = (1:th_ind);
eig_values = diag(D).^2;

%Plot principle components with index pc1 and pc2, according to the printed
%eigenvalues
%(1) If there are some nonlinear curves in the plot in any pair of index, 
%then there will be some nonlinear manifold inside the data
%(2) If data are randomly scattering in all pairs of principle components,
%then data are pretty linearly in those principle components.
PCA_score = V*D; %or X*U
title_text = [title_text,' ',num2str(pc1),' and ',num2str(pc2),' PC'];
hold on;
scatter_label2d(PCA_score(:,[pc1,pc2]),title_text,dd);
xlabel(num2str(pc1));
ylabel(num2str(pc2));
end