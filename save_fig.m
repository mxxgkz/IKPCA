%%Save figures to .jpg and .fig format

cwd = '/Users/kungangzhang/Documents/OneDrive/Northwestern/Study/Courses/Independent Study/20160720-Solve_prob_in_alg/figures/';

fig_ind = (11:16);

k0 = 0;
k = k0;
count = 0;
for i = 1:size(fig_ind,2)
    count = count + 1;
    figure(fig_ind(i));
    saveas(gca,[cwd,[num2str(k),'-',num2str(mod(i-1,3)+1)]],'jpg');
    saveas(gca,[cwd,[num2str(k),'-',num2str(mod(i-1,3)+1)]],'fig');
    if mod(count,3) ==0
        k = k + 1;
        count = 0;
    end
end