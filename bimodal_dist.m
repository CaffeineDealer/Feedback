function [bm] = bimodal_dist(x,m1,m2,sigma1,sigma2,a,b)






bm = ((a*1/2*sqrt(2*pi*(sigma1^2))) * (exp(-(x - m1).^2 / (2 * (sigma1^2))))) + ...
    ((b*1/2*sqrt(2*pi*(sigma2^2))) * (exp(-(x - m2).^2 / (2 * (sigma2^2)))));