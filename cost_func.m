function [f,g] = cost_func(Z,A,B,sigma)
%% The cost function of estimating Z given A, a, and observation X
% Calculate function value and gradient
% Input:
%   Z: Input vector. Reshape the initial N*p matrix Z into a row vector;
%   A: Coefficient matrix;
%   B: X'-a1_{1*N};
%   sigma: The sigma of Gaussian kernel
% Output:
%   f: Function value;
%   g: Gradient as a column vector
%%

% Calculate function value
[n,N] =size(A);
p = length(Z)/N;
f = 0;
Z_m = reshape(Z,[p,N])'; % Reshape the row vector Z into a N-by-p matrix
K = Gaussian_Kernel(Z_m,sigma);
for i = 1:n
    f=f+norm(B(i,:)-A(i,:)*K);
end

if nargout > 1 % graident required
    g = zeros(p*N,1);
    for k = 1:N
        for i = 1:n
            g((k-1)*p+1:k*p) = g((k-1)*p+1:k*p)+...
                             2\sigma^2*(A(i,k)*(B(i,:)-A(i,:)*K)*(K(:,k)*Z_m(k,:)-...
                             diag(K(:,k))*Z_m)+(B(i,k)-A(i,:)*K(:,k))*A(i,:)*(K(:,k)*Z_m(k,:)-...
                             diag(K(:,k))*Z_m))';
        end
    end
end
end