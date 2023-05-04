function [Rmean ctl E I] = fr2dev(Fr,meanR,theta,phi,grid)



switch meanR
    case 'end'
        FR = Fr(end,:); % Response at the end of Evolution
    case 'all'
        FR = mean(Fr);
end
rE = FR(1:grid);
rI = FR(grid+1:end);

r = []; 
r(1,:,:) = rE;
r(2,:,:) = rI;
r = squeeze(mean(r));
ctl = 0:1:180;
R = cell(1,size(ctl,2));
hx_dir = ones(1,size(theta,2)) * theta(phi);
EI_dev = round(rad2deg(abs(angdiff(hx_dir,theta))));
for j = 1:size(ctl,2)
    R{j} = [R{j} r(EI_dev == ctl(j))];
end
Rmean = cellfun(@nanmean,R);

RE = cell(1,size(ctl,2));
RI = RE;
for j = 1:size(ctl,2)
    RE{j} = [RE{j} rE(EI_dev == ctl(j))];
    RI{j} = [RI{j} rI(EI_dev == ctl(j))];
end
E = cellfun(@mean,RE);
I = cellfun(@mean,RI);


