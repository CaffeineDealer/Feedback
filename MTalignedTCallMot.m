function [Nt St errN errS dir] = MTalignedTCallMot(TCNt,TCSt,TCNs,TCSs,site)





dir = (-180:45:135);
c = 1;
for i = 1:size(TCNt,2)
    for j = 1:size(TCNt(i).MST,2)
        Nt(:,c) = TCNt(i).MST(j).MTmean;
        St(:,c) = TCSt(i).MST(j).MTmean;
        
        Ns(:,c) = TCNs(i).MST(j).MTmean;
        Ss(:,c) = TCSs(i).MST(j).MTmean;
        c = c + 1;
    end
end

ctl = [5 0;4 6;3 7;2 8;1 0];
nt = []; st = [];
ns = []; ss = [];
for j = 1:size(ctl,1)
    nt(1,:) = Nt(ctl(j,1),:);
    st(1,:) = St(ctl(j,1),:);
    
    ns(1,:) = Ns(ctl(j,1),:);
    ss(1,:) = Ss(ctl(j,1),:);
    if ~(ctl(j,2) == 0)
        nt(2,:) = Nt(ctl(j,2),:);
        st(2,:) = St(ctl(j,2),:);
        
        ns(2,:) = Ns(ctl(j,2),:);
        ss(2,:) = Ss(ctl(j,2),:);
    end
    NNt(j,:) = mean(nt,1);
    SSt(j,:) = mean(st,1);
    nt = []; st = [];
    
    NNs(j,:) = mean(ns,1);
    SSs(j,:) = mean(ss,1);
    ns = []; ss = [];
end
Nt = []; St = [];
Nt = NNt; St = SSt;

Ns = []; Ss = [];
Ns = NNs; Ss = SSs;

N = [Nt Ns];
S = [St Ss];

mN = mean(N,2)
mS = mean(S,2)
std(N,[],2)
std(S,[],2)

errN = (std(N,[],2)) / (sqrt(size(N,2)));
errS = (std(S,[],2)) / (sqrt(size(S,2)));

% ts = tinv([0.025  0.975],size(N,2)-1);
% errN = errN * ts;
% errS = errS * ts;

dir = (0:45:180);
figure
col(1) = 'k';
col(2) = 'b';

plot(dir,mN,col(1),'LineWidth',2)
hold on
plot(dir,mS,'Color',[0 .5 0],'LineWidth',2)
legend({'Control','MS'},'box','off')



shade(dir,mN+errN,sprintf('--%s',col(1)),dir,mN-errN,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,mS+errS,sprintf('--%s',[0 0.5 0]),dir,mS-errS,sprintf('--%s',[0 0.5 0]),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',[0 0.5 0])

figure
errorbar(dir,mN,errN,'Color',col(1),'LineWidth',2,'Marker','^','MarkerSize',7,'MarkerEdgeColor',col(1),'MarkerFaceColor',col(1))
hold on
errorbar(dir,mS,errS,'Color',[0 0.5 0],'LineWidth',2,'Marker','d','MarkerSize',7,'MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0 0.5 0])


% ylim([0 1])
xlim([-10 190])
% title('Mean MT TC') % MST^{Lateral} or MT
% site = 'MST^{MS site}'; % MST^{MS site} or MT
% xlabel(sprintf('Deviation from %s PD^o',site),'FontSize',15,'FontName','High Tower Text')
xlabel(['\Delta' sprintf('Pd_{%s}^{o}',site)],'FontSize',15,'FontName','High Tower Text')
ylabel('Firing Rate_{p.u.}','FontSize',15,'FontName','High Tower Text')
box off
set(gca,'XTick',[0 45 90 135 180],'XTickLabel',[0 45 90 135 180])
set(gcf,'color','w')