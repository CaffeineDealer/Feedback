function ri = vm(bw,pd,fr,spnt,dir)
%%%% Input %%%%
% params(1): Tuning Bandwidth
% params(2): Preferred Direction(PD)
% params(3): Firing Rate @ PD
% params(4): Spontaneous Firing Rate
% dir: Directions  
%%%% Output %%%%
% ri: Tuning Curve

% Written by Yavar Korkian


ri = fr .* exp(bw .* cos(dir - pd)) + spnt;
