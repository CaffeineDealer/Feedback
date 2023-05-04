function [dist] = findintsecPos(PP,ref)


dist(:,1) = pos2dist(PP(1));
dist(:,2) = pos2dist(PP(2));

for i = 1:size(dist,1)
    dist(i,3) = find(ref(:,1) == dist(i,1) & ref(:,2) == dist(i,2));    
end

