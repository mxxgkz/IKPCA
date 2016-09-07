%%Find how Gaussian kernel changes with the dimension of variation sources
%%under given sample size:
%%Z: N-by-p dimensional matrix with each row as a sample of varation
%%sources;
%%Call function Gaussian_Kernel(Z,sigma_kernel) to get the kernel matrix;
rng(2);

N = 5000;
p = [1,5,25,125];
sigma_kernel = 0.3;
pct = 0.99;

figure(1);
hold on;
for i = 1:length(p)
    Z=rand(N,p(i));
    K=Gaussian_Kernel(Z,sigma_kernel);
    e_value = eig(K);
    
    ind_thre = 1;
    sum_e_value = e_value(N);
    tot_e_value = sum(e_value);
    for j = 1:(N-1)
        if sum_e_value/tot_e_value >= pct
            break
        end
        ind_thre = ind_thre + 1;
        sum_e_value = sum_e_value + e_value(N-j);
    end
    display(['p=',num2str(p(i))]);
    display(ind_thre);
    display(e_value(N-ind_thre+1:N));
    plot((ind_thre:-1:1),e_value(N-ind_thre+1:N))
end