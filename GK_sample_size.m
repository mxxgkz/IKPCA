%%Find how Gaussian kernel changes with the sample size under a given
%%dimension of variation sources:
%%Z: N-by-p dimensional matrix with each row as a sample of varation
%%sources;
%%Call function Gaussian_Kernel(Z,sigma_kernel) to get the kernel matrix;
rng(2);

N = [500,1000,5000,10000];
p = 5;
sigma_kernel = 0.3;
pct = 0.99;

figure(1);
hold on;
for i = 1:length(N)
    Z=rand(N(i),p);
    K=Gaussian_Kernel(Z,sigma_kernel);
    e_value = eig(K);
    
    ind_thre = 1;
    sum_e_value = e_value(N(i));
    tot_e_value = sum(e_value);
    for j = 1:(N(i)-1)
        if sum_e_value/tot_e_value >= pct
            break
        end
        ind_thre = ind_thre + 1;
        sum_e_value = sum_e_value + e_value(N(i)-j);
    end
    display(['N=',num2str(N(i))]);
    display(ind_thre);
    display(e_value((N(i)-ind_thre+1):N(i)));
    plot((ind_thre:-1:1),e_value((N(i)-ind_thre+1):N(i)))
end