function th_ind = sing_th_ind(D_v,pct)
%%Find the threshold index of the singular value that corresponding to the
%%percentage of pct of variation (eigen-values)
% Input:
%   D_v: the vector of singular values
%   pct: percentage of total variation (eigen-values = square of singular
%       value
% Output:
%   th_ind: the threshold of index for the pct of variation
    sum_eig = sum(D_v.^2); %Sum of all eigen-values in PCA
    th_ind = 1; %The threshold index of principle components
    tem_sum_eig = D_v(th_ind)^2;
    n_D = length(D_v);
    while th_ind <= n_D && tem_sum_eig/sum_eig < pct
        th_ind = th_ind + 1;
        tem_sum_eig = tem_sum_eig + D_v(th_ind)^2;
    end
end