function dr = isn4(t,r,W,h,tau,k,c,n,ms)

if c == 0
    dr = (-r + k .* (W*r + h).^n) ./ tau; % Linear model equation, nonlinear if n>1
else
    dr = (-r + k .* (W*r + c.*h + ms).^n) ./ tau; % Linear model equation, nonlinear if n>1
end
% Constraint
dr(r<0) = 0;