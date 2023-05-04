function [cor] = masterOFpuppets(MT,MTms)




x = [];
for i = 1:size(MT,2)
    for j = 1:size(MT(i).df,2)
        % Find max MT pos & dir
        FrMT = MT(i).xr(j).firing;
        [x(i).Fr(j,1) I] = max(FrMT(:));
        [x(i).PDr(j,1) x(i).PDc(j,1)] = ind2sub(size(FrMT),I);
        
        x(i).MTfr(j,:) = MT(i).xr(j).firing(x(i).PDr(j,1),:);
        x(i).MTmsfr(j,:) = MTms(i).xr(j).firing(x(i).PDr(j,1),:);
    end
end

ctl = allcomb(0:1:5,0:45:180);
NMS = cell(size(ctl,1),1);
MS = NMS;
for i = 1:size(x,2)
    for j = 1:size(MT(i).df,2)
        out = zeros(8,9,3);
        outms = out;
        out(:,:,1) = MT(i).xr(j).firing;
        outms(:,:,1) = MTms(i).xr(j).firing;        
        % Dir
        pd = x(i).PDr(j,1);
        d = findintsec(indx2dir(pd),ctl(1:5,2));
        ddist(:,1) = dir2indx(d(:,1));
        ddist(:,2) = dir2indx(d(:,2));
        ddist(:,3) = ctl(1:5,2);
        for z = 1:size(ddist,1)
            out(unique(ddist(z,1:2)),:,2) = ddist(z,3);
        end
        % Pos
        pp = x(i).PDc(j,1);
        pdist = pos2dist(pp);
        for w = 1:9
            out(:,w,3) = pdist(w,1);
        end
        outms(:,:,2:3) = out(:,:,2:3); 
        %
        for u = 1:size(ctl,1)
            val = out(out(:,:,2) == ctl(u,2) & out(:,:,3) == ctl(u,1));
            valms = outms(out(:,:,2) == ctl(u,2) & out(:,:,3) == ctl(u,1));
            NMS{u,1} = [NMS{u,1} mean(val)];
            MS{u,1} = [MS{u,1} mean(valms)];
            val = [];
        end
    end    
end

NMS = cellfun(@nanmean,NMS);
MS = cellfun(@nanmean,MS);
firing = MS - NMS;
cor = flip(reshape(firing,5,6)');


