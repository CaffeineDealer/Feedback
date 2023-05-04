function dr = isn2E(t,r,W,h,tau,c,k,n,a,ms)


dr = (-r + k .* (W*r + c.*h + a.*ms).^n) ./ tau; % Linear model equation, nonlinear if n>1

% Constraint
dr(r<0) = 0;