function z = Gauss2D(x,y,x0,y0,r,plt)
%%%% Input %%%%
% x: 
% y: 
% x0: center of Gaussian
% y0: center of Gaussian
% r:  width 
%%%% Output %%%%
% z: gauss

% Written by Yavar Korkian
% Last Date Modified: 19.09.2021




z = exp((-((x - x0).^2 + (y - y0).^2)) / (2 * (r)^2));

if plt == 1
%     figure
%     mesh(x,y,z)
%     colormap([0 0 0])
%     figure
    v = 0:0.01:1;
    contour(x,y,z,v)
    axis square
end