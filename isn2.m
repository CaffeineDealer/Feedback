function dr = isn2(t,r,W,h,tau,c,k,n)


dr = (-r + k .* (W*r + c.*h).^n) ./ tau; % Linear model equation, nonlinear if n>1

% Constraint
% dr(r<0) = 0;