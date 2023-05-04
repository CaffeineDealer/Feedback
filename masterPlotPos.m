function [cor] = masterPlotPos(MST,MT,MTms,type)





x = [];
for i = 1:size(MT,2)
    for z = 1:size(MST(i).df,2)
        FrMST = MST(i).xr(z).firing;
        [frr I] = max(FrMST(:)); % x(i).inf(z).Fr(j,1)
        for j = 1:size(MT(i).df,2)
            % Find max MST pos & dir
            [x(i).inf(z).PDr(j,1) x(i).inf(z).PDc(j,1)] = ind2sub(size(FrMST),I);
            % Find max MT pos & dir
            FrMT = MT(i).xr(j).firing;
            [x(i).inf(z).Fr(j,2) II] = max(FrMT(:));
            [x(i).inf(z).PDr(j,2) x(i).inf(z).PDc(j,2)] = ind2sub(size(FrMT),II);
            % 
            x(i).MTfr(z,j,:) = MT(i).xr(j).firing(x(i).inf(z).PDr(j,2),:); 
            x(i).MTmsfr(z,j,:) = MTms(i).xr(j).firing(x(i).inf(z).PDr(j,2),:); 
        end
    end
end

ctl = allcomb(0:1:5,0:1:5);
Fr = [];
NMS = cell(size(ctl,1),1);
MS = NMS;
for i = 1:size(x,2)
    for z = 1:size(MST(i).df,2)
        for j = 1:size(MT(i).df,2)
            Frmt = squeeze(x(i).MTfr(z,j,:));
            Frmtms = squeeze(x(i).MTmsfr(z,j,:));
            
            mstPP = x(i).inf(z).PDc(j,1);
            mtPP = x(i).inf(z).PDc(j,2);
            
            dist = findintsecPos([mstPP,mtPP],ctl);
            upos = unique(dist(:,3));
            for w = 1:size(upos,1)
                NMS{upos(w)} = [NMS{upos(w)} mean(Frmt(find(dist(:,3) == upos(w))))]; 
                MS{upos(w)} = [MS{upos(w)} mean(Frmtms(find(dist(:,3) == upos(w))))]; 
            end
        end
    end
end

nms = cellfun(@mean,NMS);
ms = cellfun(@mean,MS);
firing = ms - nms;
% firing = ms;
cor = reshape(firing,6,6)';

figure
% colormap(summer)
imagesc(cor)
set(gca,'YDir','Normal')
caxis([min(cor(:)) abs(min(cor(:)))])
% colormap jet
title(sprintf('%s',type))
xlabel 'Stim Pos from MT PP' %MST^{Lateral}
ylabel 'Stim Pos from MST PP' %^{MS Site}
xl = 0:1:5;
yl = xl;
set(gca,'XTick',[1:1:6],'XTickLabel',xl)
set(gca,'YTick',[1:1:6],'YTickLabel',yl)
colorbar
