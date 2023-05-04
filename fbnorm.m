function [R] = fbnorm(Ii,Ij,params)





% Written by Yavar Korkian
% Date: 02.09.2021
% Last Update: 02.09.2021


if size(Ij,1) == 1
    error('Yo! I am expecting a vector of firing rate')
end

sigma = params.sigma;
w = params.w;
M = params.M;


n = Ii;
d = sum((w .* Ij),1);

R = (M * n) ./ (d + sigma);
