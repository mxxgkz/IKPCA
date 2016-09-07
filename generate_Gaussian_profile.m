function [X] = generate_Gaussian_profile(Z,N,sigma_data,l)
%%Generate data as profiles of Gaussian distribution with given sigma_data,
%%but different center location
%%Z: N*2 matrix with each row as a 2-D center location of the Gaussian
%%profile
%%N: the sample size of those generated "images"
%%sigma_data: the sigma of those Gaussian profiles
%%l: the largest index of pixel in one side of images of those Gaussian
%%profiles (index starts from zero)

rng(2);

X = zeros(N,(l+1)^2);
center_loc = Z;

for i = 1:N
    for j = 0:l
        for k = 0:l
            X(i,j*(l+1)+k+1) = exp(-((j-center_loc(i,1))^2+(k-center_loc(i,2))^2)/2/sigma_data^2);
        end
    end
end

end