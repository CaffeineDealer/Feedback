function [dfiring X] = fr_rank(datapath,bin,MS)





inf = load([datapath sprintf('N%s.mat','trans')]);
MTnmsT = inf.info;
clear inf

inf = load([datapath sprintf('N%s.mat','spiral')]);
MTnmsS = inf.info;
clear inf

inf = load([datapath sprintf('S%s.mat','trans')]);
MTmsT = inf.info;
clear inf

inf = load([datapath sprintf('S%s.mat','spiral')]);
MTmsS = inf.info;
clear inf

if MS == 1
    MTnmsT(7) = [];
    MTnmsT(4) = [];
    MTnmsT(1) = [];
    
    MTnmsS(7) = [];
    MTnmsS(4) = [];
    MTnmsS(1) = [];
    
    
    MTmsT(7) = [];
    MTmsT(4) = [];
    MTmsT(1) = [];
    
    MTmsS(7) = [];
    MTmsS(4) = [];
    MTmsS(1) = [];
end

c = 1;
for i = 1:size(MTnmsT,2)
    for j = 1:size(MTnmsT(i).xr,2)
        a(:,1) = reshape(vertcat(MTnmsT(i).xr(j).firing,MTnmsS(i).xr(j).firing),144,1);
        a(:,2) = reshape(vertcat(MTmsT(i).xr(j).firing,MTmsS(i).xr(j).firing),144,1);
        X(:,:,c) = a;
        a = flip(sortrows(a,1));
        nms(:,c) = a(:,1);
        ms(:,c) = a(:,2);
        c = c + 1;
        a = [];
    end
end

fr = [];
fr = ms - nms;
dfiring = flip(mean(fr,2));
gmin = abs(min(dfiring));
gmax = abs(max(dfiring));
glb = max(gmin,gmax);
gmax = 144;
offset = 0.2;
h = figure;
patch([0 gmax+offset gmax+offset 0],[0 0 glb+offset glb+offset],[0.5 0.5 0.5],'FaceAlpha',.2)
patch([0 gmax+offset gmax+offset 0],[0 0 -glb-offset -glb-offset],[0.5 0.5 0.5],'FaceAlpha',.4)
hold on
bar(dfiring,.2,'k')
set(gca,'XTick',[])
xlim([1-offset gmax+offset])
ylim([-glb-offset glb+offset])
xlabel('Worst \rightarrow Best','FontSize',15,'FontName','High Tower Text')
ylabel('\DeltaFr_{p.u.}','FontSize',15,'FontName','High Tower Text')

% ctl = linspace(0,1,bin)';
% for i = 1:size(nms,2)
%     for j = 1:size(ctl,1)-1
%         if j ~= size(ctl,1)-1
%             fr(i,j,1) = mean(nms(nms(:,i) >= ctl(j) & nms(:,i) < ctl(j+1),i));
%             fr(i,j,2) = mean(ms(nms(:,i) >= ctl(j) & nms(:,i) < ctl(j+1),i));
%         elseif j == size(ctl,1)-1
%             fr(i,j,1) = mean(nms(nms(:,i) >= ctl(j) & nms(:,i) <= ctl(j+1),i));
%             fr(i,j,2) = mean(ms(nms(:,i) >= ctl(j) & nms(:,i) <= ctl(j+1),i));
%         end
%     end    
% end

% frErr = squeeze(fr(:,:,2) - fr(:,:,1));
% err = (nanstd(frErr)) / (sqrt(size(frErr,1)));
% 
% FR = (squeeze(nanmean(fr,1)));
% firing = FR(:,2) - FR(:,1);
% gmin = abs(min(firing));
% gmax = abs(max(firing));
% glb = max(gmin,gmax);
% gmax = size(ctl,1) - 1;
% offset = 0.2;
% %%%%%%%% 
% figure
% patch([0 gmax+offset gmax+offset 0],[0 0 glb+offset glb+offset],[0.5 0.5 0.5],'FaceAlpha',.2)
% patch([0 gmax+offset gmax+offset 0],[0 0 -glb-offset -glb-offset],[0.5 0.5 0.5],'FaceAlpha',.4)
% hold on
% bar(firing,0.2,'k')
% errorbar(firing,err,'Color',[0, 0.5, 0],'LineWidth',2,'LineStyle','none')
% xlim([1-offset gmax+offset])
% ylim([-glb-offset glb+offset])
% xt = floor((ctl(2:end,:)*100));
% xl = (1:size(ctl,1)-1)';
% % set(gca,'XTick',xl,'XTickLabel',xt)
% set(gca,'XTick',[])
% % title 'Firing Rate Breakdown'
% % xlabel 'Top %'
% % ylabel 'Norm {\Delta}FR = Fr_{MT}^{MS} - Fr_{MT}^{NMS}'
% xlabel('Worst \rightarrow Best','FontSize',15,'FontName','High Tower Text')
% ylabel('\DeltaFr_{p.u.}','FontSize',15,'FontName','High Tower Text')