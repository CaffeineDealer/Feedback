% PD difference
clear all; clc
datapath = 'D:\MT_MST\Microstim\PSTH\MUA_SUA\Norm\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
[pdDiff Fr] = prefDiff(datapath,20);
%%
clear all; clc
dp = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
Nt = load([dp 'Ntrans.mat']);
Ns = load([dp 'Nspiral.mat']);
St = load([dp 'Strans.mat']);
Ss = load([dp 'Sspiral.mat']);


Nt.info(7) = []; Nt.info(4) = []; Nt.info(1) = [];
Ns.info(7) = []; Ns.info(4) = []; Ns.info(1) = [];

St.info(7) = []; St.info(4) = []; St.info(1) = [];
Ss.info(7) = []; Ss.info(4) = []; Ss.info(1) = [];

N1 = struct2cell(Nt);
N2 = struct2cell(Ns);
Nt = [N1{1,1}(:).xr];
Ns = [N2{1,1}(:).xr];

S1 = struct2cell(St);
S2 = struct2cell(Ss);
St = [S1{1,1}(:).xr];
Ss = [S2{1,1}(:).xr];

[Vxt Vyt diffVt] = RFupdated(St,Nt);
[Vxs Vys diffVs] = RFupdated(Ss,Ns);

for i = 1:size(Vxt,1)
    dt(i,1) = abs(Edist(Vxt(i,:),Vyt(i,:))); % Translation
    dt(i,2) = abs(Edist(Vxs(i,:),Vys(i,:))); % Spirals
end
%% MT Tuning curves aligned with MST: Part II.II
clear all; clc
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
site = 'MT';
MS = 0; % 0 if MT, 1 if MST

switch site
    case 'MST'
        TCNt = PDFalign(datapath,'N','trans',0,MS);
        TCSt = PDFalign(datapath,'S','trans',0,MS);
        
        TCNs = PDFalign(datapath,'N','spiral',0,MS);
        TCSs = PDFalign(datapath,'S','spiral',0,MS);
    case 'MT'        
        [TCNt Xnt] = PDFalignMT(datapath,'N','trans',0,MS);
        [TCSt Xst] = PDFalignMT(datapath,'S','trans',0,MS);
        
        [TCNs Xns] = PDFalignMT(datapath,'N','spiral',0,MS);
        [TCSs Xss] = PDFalignMT(datapath,'S','spiral',0,MS);
end

[N S errN errS dirt] = MTalignedTCallMot(TCNt,TCSt,TCNs,TCSs,site);

% for i = 1:size(N,1)
%     [h(i),p(i)] = ttest(N(i,:),S(i,:));
%     if p(i) <= 0.05
%         stat{i} = 'Significant';
%     else
%         stat{i} = 'Non-Significant';
%     end
%     sprintf('%s=%.3f',stat{i},p(i))
% end
n = reshape(N,size(N,1)*size(N,2),1);
s = reshape(S,size(S,1)*size(S,2),1);
[hns pns] = ttest(n,s)
MeanN = mean(n)
stdN = std(n)
MeanS = mean(s)
stdS = std(s)
%
Xnms = [Xnt Xns];
Xms = [Xst Xss];
for i = 1:5
    pwpNMS(:,i) = pwproduct(Xnms(i,:));
    pwpMS(:,i) = pwproduct(Xms(i,:));
end
pwpNMS = reshape(pwpNMS,size(pwpNMS,1)*size(pwpNMS,2),1);
pwpMS = reshape(pwpMS,size(pwpMS,1)*size(pwpMS,2),1);
pwp = pwpMS - pwpNMS;

h = figure;
e = 0.05;
lb = min(pwp)- e;
ub = max(pwp) + e;
histogram(pwp,'FaceColor','[0.75, 0.75, 0]','EdgeColor','[0.75, 0.75, 0]','Normalization','count','FaceAlpha',.8)
set(gca,'YScale','log')
xlabel('\DeltaPairwise Product','FontSize',15,'FontName','High Tower Text')
box 'off'
h.Color = 'w';
xlim([lb ub])
[h,p] = ttest(pwp);
[bf,bc] = bimodalitycoeff(pwp);
%% Position @ PDirection
clear all; clc
% datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
MS = 0;
[MTt NMSmtT MSmtT bl] = PosAtPD(datapath,'trans',MS,'MT');
[MSTt NMSmstT MSmstT bl] = PosAtPDMST(datapath,'trans',MS,'MST');

[MTs NMSmtS MSmtS bl] = PosAtPD(datapath,'spiral',MS,'MT');
[MSTs NMSmstS MSmstS bl] = PosAtPDMST(datapath,'spiral',MS,'MST');



for i = 1:size(NMSmtT,2)
    a = NMSmtT(~isnan(NMSmtT(:,i)),i)';
    b = NMSmtS(~isnan(NMSmtS(:,i)),i)';
    n = [a b]';
    
    c = MSmtT(~isnan(MSmtT(:,i)),i)';
    d = MSmtS(~isnan(MSmtS(:,i)),i)';
    s = [c d]';
    
    errNmt(1,i) = (std(n)) / (sqrt(size(n,1))); 
    errSmt(1,i) = (std(s)) / (sqrt(size(s,1))); 
    n = []; s = [];
    a = []; b = []; c = []; d = [];
end

NMS = [NMSmtT NMSmtS];
NMS = reshape(NMS,size(NMSmtT,1),size(NMSmtT,2),2);
NMS = squeeze(nanmean(NMS,3));
NMS = squeeze(nanmean(NMS));

MS = [MSmtT MSmtS];
MS = reshape(MS,size(MSmtT,1),size(MSmtT,2),2);
MS = squeeze(nanmean(MS,3));
MS = squeeze(nanmean(MS));

for i = 1:size(NMSmstT,2)
    a = NMSmstT(~isnan(NMSmstT(:,i)),i)';
    b = NMSmstS(~isnan(NMSmstS(:,i)),i)';
    n = [a b]';
    
    c = MSmstT(~isnan(MSmstT(:,i)),i)';
    d = MSmstS(~isnan(MSmstS(:,i)),i)';
    s = [c d]';
    
    errNmst(1,i) = (std(n)) / (sqrt(size(n,1))); 
    errSmst(1,i) = (std(s)) / (sqrt(size(s,1))); 
    n = []; s = [];
    a = []; b = []; c = []; d = [];
end

NMSmst = [NMSmstT NMSmstS];
NMSmst = reshape(NMSmst,size(NMSmstT,1),size(NMSmstT,2),2);
NMSmst = squeeze(nanmean(NMSmst,3));
NMSmst = squeeze(nanmean(NMSmst));

MSmst = [MSmstT MSmstS];
MSmst = reshape(MSmst,size(MSmstT,1),size(MSmstT,2),2);
MSmst = squeeze(nanmean(MSmst,3));
MSmst = squeeze(nanmean(MSmst));

figure
site = 'MST^{Lateral}';
dir = 1:6;
col(1) = 'k';
col(2) = 'b';
plot(dir,NMS,'Color',col(1),'LineWidth',2); hold on
plot(dir,MS,'Color',col(2),'LineWidth',2)
title 'Position'
xlabel(sprintf('Position Deviation from %s',site))
ylabel 'Normalized Firing Rate'
legend({'NMS','MS'})
xticks([1 2 3 4 5 6])
xticklabels({'0','1','2','3','4','5'})

shade(dir,NMS+errNmt,sprintf('--%s',col(1)),dir,NMS-errNmt,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,MS+errSmt,sprintf(':%s',col(2)),dir,MS-errSmt,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

line([min(dir) max(dir)],[bl bl],'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--')

figure
site = 'MST^{MS site}';
dir = 1:6;
col(1) = 'k';
col(2) = 'b';
plot(dir,NMSmst,'Color',col(1),'LineWidth',2); hold on
plot(dir,MSmst,'Color',col(2),'LineWidth',2)
title 'Position'
xlabel(sprintf('Position Deviation from %s',site))
ylabel 'Normalized Firing Rate'
legend({'NMS','MS'})
xticks([1 2 3 4 5 6])
xticklabels({'0','1','2','3','4','5'})

shade(dir,NMSmst+errNmst,sprintf('--%s',col(1)),dir,NMSmst-errNmst,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,MSmst+errSmst,sprintf(':%s',col(2)),dir,MSmst-errSmst,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))

line([min(dir) max(dir)],[bl bl],'LineWidth',2,'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--')


%% Master Plot(Fr & PD): MT vs MST 
clear all; clc
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
type = 'trans';
ty = 'spiral';
MS = 1;
inf = load([datapath sprintf('N%sMS.mat',type)]); 
MSTt = inf.info; 
clear info inf
inf = load([datapath sprintf('N%s.mat',type)]); 
MTtnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',type)]); 
MTtms = inf.info; 
clear info inf
%
inf = load([datapath sprintf('N%sMS.mat',ty)]); 
MSTs = inf.info; 
clear info inf
inf = load([datapath sprintf('N%s.mat',ty)]); 
MTsnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',ty)]); 
MTsms = inf.info; 
clear info inf

for i = 1:size(MTtnms,2)
    for j = 1:size(MTtnms(i).xr,2)
        MTtms(i).df(:,j) = MTtms(i).xr(j).firing(:,MTtnms(i).pt(j));
        MTsms(i).df(:,j) = MTsms(i).xr(j).firing(:,MTsnms(i).pt(j));
    end
end

if MS == 1
    MTtnms(7) = [];
    MTtnms(4) = [];
    MTtnms(1) = [];
    MTtms(7) = [];
    MTtms(4) = [];
    MTtms(1) = [];
    MTsnms(7) = [];
    MTsnms(4) = [];
    MTsnms(1) = [];
    MTsms(7) = [];
    MTsms(4) = [];
    MTsms(1) = [];
end

[cort] = masterPlot(MSTt,MTtnms,MTtms,'Translation'); %x
[cors] = masterPlot(MSTs,MTsnms,MTsms,'Spiral');
cor = [cort cors];
cor = reshape(cor,size(cort,1),size(cors,1),2);
cor = squeeze(mean(cor,3));

[cortP] = masterPlotPos(MSTt,MTtnms,MTtms,'Translation');
[corsP] = masterPlotPos(MSTs,MTsnms,MTsms,'Spiral');
corP = [cortP corsP];
corP = reshape(corP,size(cortP,1),size(corsP,1),2);
corP = squeeze(mean(corP,3));

% Plot Master Direction & Master Position
h1 = figure;
colormap gray
imagesc(cor)
caxis([min(cor(:)) abs(min(cor(:)))])
box off
h1.Color = 'w';
xl = 0:45:180;
yl = 180:-45:0;
% title('Direction','FontSize',15,'FontName','High Tower Text')
xlabel('Stim Dir from MT PD^{o}','FontSize',15,'FontName','High Tower Text') %MST^{Lateral}
ylabel('Stim Dir from MST PD^{o}','FontSize',15,'FontName','High Tower Text') %^{MS Site}
set(gca,'XTick',[1:1:5],'XTickLabel',xl)
set(gca,'YTick',[1:1:5],'YTickLabel',yl)
colorbar


h2 = figure;
colormap pink
imagesc(corP)
caxis([min(corP(:)) abs(min(corP(:)))])
box off
h2.Color = 'w';
xl = 0:1:5;
yl = xl;
% title('Position','FontSize',15,'FontName','High Tower Text')
xlabel('Stim Pos from MT PP','FontSize',15,'FontName','High Tower Text') %MST^{Lateral}
ylabel('Stim Pos from MST PP','FontSize',15,'FontName','High Tower Text') %^{MS Site}
set(gca,'XTick',[1:1:6],'XTickLabel',xl)
set(gca,'YTick',[1:1:6],'YTickLabel',yl)
set(gca,'YDir','Normal')
colorbar
%% Sooo Master
clear all; clc
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
type = 'trans';
ty = 'spiral';
MS = 0;
% Translation
inf = load([datapath sprintf('N%s.mat',type)]); 
MTtnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',type)]); 
MTtms = inf.info; 
clear info inf
% Spirals
inf = load([datapath sprintf('N%s.mat',ty)]); 
MTsnms = inf.info; 
clear info inf
inf = load([datapath sprintf('S%s.mat',ty)]); 
MTsms = inf.info; 
clear info inf

for i = 1:size(MTtnms,2)
    for j = 1:size(MTtnms(i).xr,2)
        MTtms(i).df(:,j) = MTtms(i).xr(j).firing(:,MTtnms(i).pt(j));
        MTsms(i).df(:,j) = MTsms(i).xr(j).firing(:,MTsnms(i).pt(j));
    end
end

if MS == 1
    MTtnms(7) = [];
    MTtnms(4) = [];
    MTtnms(1) = [];
    MTtms(7) = [];
    MTtms(4) = [];
    MTtms(1) = [];
    MTsnms(7) = [];
    MTsnms(4) = [];
    MTsnms(1) = [];
    MTsms(7) = [];
    MTsms(4) = [];
    MTsms(1) = [];
end

corT = masterOFpuppets(MTtnms,MTtms);
corS = masterOFpuppets(MTsnms,MTsms);

cor = [corT corS];
cor = reshape(cor,size(corT,1),size(corT,2),2);
cor = squeeze(mean(cor,3));

h3 = figure;
colormap bone
imagesc(cor)
caxis([min(cor(:)) abs(min(cor(:)))])
box off
h3.Color = 'w';
% title 'Soooo Master'
xlabel('Stim Dir from MT PD^{o}','FontSize',15,'FontName','High Tower Text') %MST^{Lateral}
ylabel('Stim Pos from MT PP','FontSize',15,'FontName','High Tower Text') %^{MS Site}
xl = 0:45:180;
yl = 5:-1:0;
set(gca,'XTick',[1:1:5],'XTickLabel',xl)
set(gca,'YTick',[1:1:6],'YTickLabel',yl)
colorbar
%% Firing Rate Rank Breakdown
clear all; clc
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'D:\MT_MST\Microstim\PSTH\ms_rmvd_170toEnd\Normalized\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
MS = 0;
bin = 144;
[firing X] = fr_rank(datapath,bin,MS);
%
x1 = squeeze(X(:,1,:)); %NMS
x2 = squeeze(X(:,2,:)); %MS
pwp = [];
for i = 1:144
    pwp(i,:,1) = pwproduct(x1(i,:)); %NMW
    pwp(i,:,2) = pwproduct(x2(i,:)); %MS
end

pmean = squeeze(mean(pwp,1));
PWP = pmean(:,2) - pmean(:,1);
ubx = max(pmean(:,1));
uby = max(abs(PWP));
h0 = figure; 
scatter(pmean(:,1),abs(PWP),'.','k'); hold on
h1 = refline(1,0)
h0.Color = 'w';
h1.Color = 'r';
h1.LineWidth = 2;
xlim([0 0.13])
ylim([0 0.07])
xlabel('C_{ij}^{NMS}','FontSize',15,'FontName','High Tower Text')
ylabel('|\DeltaC_{ij}|','FontSize',15,'FontName','High Tower Text')
%%
clear all; clc
load('D:\MT_MST\Microstim\Model\mt.mat','TCNt');
Nt1 = TCNt;
load('D:\MT_MST\Microstim\Model\mst.mat','TCNt');
Nt2 = TCNt;

load('D:\MT_MST\Microstim\Norm\MUA_SUA\Ntrans.mat')
info(7) = [];
info(4) = [];
info(1) = [];
dMT = {info.Dir};
dMT = cat(1,dMT{:});
dMT(:,2) = []; 
info = [];
load('D:\MT_MST\Microstim\Norm\MUA_SUA\NtransMS.mat')
dMST = {info.Dir};
dMST = cat(1,dMST{:});
dMST(:,2) = []; 
clear info
load('D:\MT_MST\Microstim\Model\pd.mat')

mt1 = [];
mt2 = [];
for i = 1:size(Nt1,2)
    a = size(Nt2(i).MST,2);
    csize(i,1) = a;
    csize(i,2) = size(Nt1(i).MST.bMT,2);
    for j = 1:a
        mt1 = [mt1 Nt1(i).MST.bMT];
    end
    mt2 = [mt2 Nt2(i).MST.bMT]; 
end

ctl = [5 0;4 6;3 7;2 8;1 0];
x1 = [];
x2 = [];
for i = 1:size(mt1,2)
    for j = 1:size(ctl,1)
        x1(1,:) = mt1(ctl(j,1),i);
        x2(1,:) = mt2(ctl(j,1),i);
        if ~(ctl(j,2) == 0)
            x1(2,:) = mt1(ctl(j,2),i);
            x2(2,:) = mt2(ctl(j,2),i);
        end
        X1(j,i) = mean(x1,1);
        X2(j,i) = mean(x2,1);
        x1 = [];
        x2 = [];
    end
end

ctl = csize(:,1) .* csize(:,2);
cnt = zeros(sum(csize(:,1)),2);
a = 0;
for i = 1:size(ctl,1)
    for j = 1:csize(i,1)
        a = a + 1;
        if a == 1
            cnt(a,1) = 1;
            cnt(a,2) = csize(i,2) * j;
        else 
            cnt(a,1) = cnt(a-1,2) + 1;
            cnt(a,2) = csize(i,2) + cnt(a,1) - 1;
        end
    end    
end
%%
R = [];
for i = 1:size(cnt,1)
    ncell = cnt(i,2) - cnt(i,1) + 1;%size(X1,2);
    params.M = 1;
    r = (1 - 0.5) .* rand(ncell,1) + 0.5;
%     r = zeros(ncell,1) + 0.1;
    params.w = r; %/sum(r);
    params.sigma = 0.00005;
    
    res = fbnorm(X1(:,cnt(i,1):cnt(i,2)),X2(:,cnt(i,1):cnt(i,2)),params);
    R = [res R];
end

figure
for i = 1:size(R,2)
    subplot(14,14,i)
    plot(X1(:,i),'b','LineWidth',0.6); hold on
    plot(X2(:,i),'k','LineWidth',0.6)
    plot(R(:,i),'r','LineWidth',1.2)
end
legend({'Observed MT^{aligned with itself}','Observed MT^{aligned with MST}','Normalized MT'})

figure
plot(mean(X1,2),'b','LineWidth',1.5); hold on
plot(mean(X2,2),'k','LineWidth',1.5)
plot(mean(R,2),'r','LineWidth',2)
legend({'Observed MT^{aligned with itself}','Observed MT^{aligned with MST}','Normalized MT'})
%%
Y(:,1) = x;
Y(:,2) = R;
Y = sort(Y,1,'descend');

figure
subplot(1,2,1)
scatter(Y(:,1),Y(:,2),'b','filled')
set(gca,'YScale','log'); grid on
% xlabel 'Cell Resp before Norm'
% ylabel 'Cell Resp after Norm'
% refline(1,0)
% subplot(1,4,2)
% plot(Y(:,2),'LineWidth',2,'Color','k')
% grid on
% subplot(1,4,3)
% plot(Y(:,1),'LineWidth',2,'Color','r')
% grid on
subplot(1,2,2)
plot(Y(:,1),'LineWidth',2,'Color','k'); hold on
plot(Y(:,1)-Y(:,2),'LineWidth',2,'Color','r')
plot(Y(:,2),'LineWidth',2,'Color','b')
grid on
legend({'In','Diff','Norm In'})
%%
clear all; clc
% datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
MS = '';
inf = load([datapath sprintf('Ntrans%s.mat',MS)]); 
trans = inf.info; 
clear info inf
inf = load([datapath sprintf('Nspiral%s.mat',MS)]);
spiral = inf.info;
clear info inf

y = [];
[y(1).r y(1).p y(1).ch] = CorrDist(trans);
[y(2).r y(2).p y(2).ch] = CorrDist(spiral);

figure
c = 1;
for i = 1:size(y(c).r,1)
    plot(y(c).ch{i},y(c).r{i},'*b','MarkerSize',5); hold on 
    plot(y(c+1).ch{i},y(c+1).r{i},'*k','MarkerSize',5); hold on 
end
a = max(cellfun(@max,y(1).ch));
xlim([0 a + 0.5])
xticks(1:a)
xticklabels(150:150:a*150)
legend({'Trans','Spirals'})
ylim([-1 1])
xlabel 'Distance (\mum)'
ylabel 'Correlation Coefficient'
title 'Intersite Correlation of MST Selectivity'
%%
clear all; clc
datapath = 'D:\MT_MST\Microstim\';
load([datapath 'MST.mat'])

y = [];
[y(1).r y(1).p y(1).ch] = CorrDist(trans);
[y(2).r y(2).p y(2).ch] = CorrDist(spiral);

figure
c = 1;
for i = 1:size(y(c).r,1)
    subplot(1,2,1)
    plot(y(c).ch{i},y(c).r{i},'*b','MarkerSize',5); hold on 
end
a = max(cellfun(@max,y(c).ch));
xlim([0 a + 0.5])
xticks(1:a)
xticklabels(150:150:a*150)
legend({'Trans'})
ylim([-1 1])
xlabel 'Distance (\mum)'
ylabel 'Correlation Coefficient'


c = 2;
for i = 1:size(y(c).r,1)
    subplot(1,2,2)
    plot(y(c).ch{i},y(c).r{i},'*k','MarkerSize',5); hold on 
end
a = max(cellfun(@max,y(c).ch));
xlim([0 a + 0.5])
xticks(1:a)
xticklabels(150:150:a*150)
legend({'Spiral'})
ylim([-1 1])
xlabel 'Distance (\mum)'
suptitle 'Intersite Correlation of MST Selectivity'
%%
clear all; clc
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% datapath = 'D:\MT_MST\Microstim\PSTH\MUA-MST-V1V2\Norm\';
MS = 1;
inf = load([datapath 'NtransMS.mat']); 
mstT = inf.info; 
clear info inf
inf = load([datapath 'Ntrans.mat']);
mtT = inf.info;
clear info inf

inf = load([datapath 'NspiralMS.mat']); 
mstS = inf.info; 
clear info inf
inf = load([datapath 'Nspiral.mat']);
mtS = inf.info;
clear info inf

if MS == 1
    mtT(7) = [];
    mtT(4) = [];
    mtT(1) = [];
    mtS(7) = [];
    mtS(4) = [];
    mtS(1) = [];
end

wT = CorrDistMSTMT(mstT,mtT);
wS = CorrDistMSTMT(mstS,mtS);
save([datapath 'w.mat'],'wT','wS')
%%
clear all; clc
load('D:\MT_MST\Microstim\Model\mt.mat','TCNt');
Nt1 = TCNt;
clear TCNt;
load('D:\MT_MST\Microstim\Model\mst.mat','TCNt');
Nt2 = TCNt;
clear TCNt;
datapath = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
load([datapath 'w.mat'])

mt1 = [];
mt2 = [];
for i = 1:size(Nt1,2)
    a = size(Nt2(i).MST,2);
    csize(i,1) = a;
    csize(i,2) = size(Nt1(i).MST.bMT,2);
    for j = 1:a
        mt1 = [mt1 Nt1(i).MST.bMT];
    end
    mt2 = [mt2 Nt2(i).MST.bMT]; 
end

ctl = [5 0;4 6;3 7;2 8;1 0];
x1 = [];
x2 = [];
for i = 1:size(mt1,2)
    for j = 1:size(ctl,1)
        x1(1,:) = mt1(ctl(j,1),i);
        x2(1,:) = mt2(ctl(j,1),i);
        if ~(ctl(j,2) == 0)
            x1(2,:) = mt1(ctl(j,2),i);
            x2(2,:) = mt2(ctl(j,2),i);
        end
        X1(j,i) = mean(x1,1);
        X2(j,i) = mean(x2,1);
        x1 = [];
        x2 = [];
    end
end

ctl = csize(:,1) .* csize(:,2);
cnt = zeros(sum(csize(:,1)),2);
a = 0;
for i = 1:size(ctl,1)
    for j = 1:csize(i,1)
        a = a + 1;
        if a == 1
            cnt(a,1) = 1;
            cnt(a,2) = csize(i,2) * j;
        else 
            cnt(a,1) = cnt(a-1,2) + 1;
            cnt(a,2) = csize(i,2) + cnt(a,1) - 1;
        end
    end    
end

w = [];
for i = 1:size(wT,2)
    w = [w ; vertcat(wT(i).r(:))];
end

w = (w);
w = ones(194,1)*1;
% tit = 'w = corrcoef (MST,MT)';
tit = 'w = -1';

R = [];
for i = 1:size(cnt,1)
    ncell = cnt(i,2) - cnt(i,1) + 1;%size(X1,2);
    params.M = 1;
    r = (1 - 0.5) .* rand(ncell,1) + 0.5;
%     r = zeros(ncell,1) + 0.1;
    params.w = r; %/sum(r);
    params.w = w(cnt(i,1):cnt(i,2));
    params.sigma = 0.00005;
    
    res = fbnorm(X1(:,cnt(i,1):cnt(i,2)),X2(:,cnt(i,1):cnt(i,2)),params);
    R = [res R];
end

figure
dir = 0:45:180;
plot(dir,mean(X1,2),'b','LineWidth',1.5); hold on
plot(dir,mean(X2,2),'k','LineWidth',1.5)
plot(dir,mean(R,2),'r','LineWidth',2)
legend({'Observed MT^{aligned with itself}','Observed MT^{aligned with MST}','Normalized MT'})
xlabel 'Deviation from MT PD^{\circ}'
ylabel 'Normalized Firing Rate'
xlim([-5 185])
title(tit)
%%
figure
for i = 1:size(R,2)
    subplot(14,14,i)
    plot(X1(:,i),'b','LineWidth',0.6); hold on
    plot(X2(:,i),'k','LineWidth',0.6)
    plot(R(:,i),'r','LineWidth',1.2)
end
legend({'Observed MT^{aligned with itself}','Observed MT^{aligned with MST}','Normalized MT'})
%% PSTH Paper
clear all; clc
datapath = 'D:\MT_MST\Microstim\PSTH\ms_removed_all\';
Ntrans = load([datapath 'Ntrans.mat'],'info');
Nspiral = load([datapath 'Nspiral.mat'],'info');
Strans = load([datapath 'Strans.mat'],'info');
Sspiral = load([datapath 'Sspiral.mat'],'info');

N = [horzcat(Ntrans.info.sp) horzcat(Nspiral.info.sp)];
S = [horzcat(Strans.info.sp) horzcat(Sspiral.info.sp)];

for i = 1:size(N,2)
    sizen(i,1) = size(N(i).spktrain,1);
    sizen(i,2) = size(N(i).spktrain_bl,1);
    sizes(i,1) = size(S(i).spktrain,1);
    sizes(i,2) = size(S(i).spktrain_bl,1);
end
spT = min([min(sizen(:,1)) min(sizes(:,1))]);
spB = min([min(sizen(:,2)) min(sizes(:,2))]);

for i = 1:size(N,2)
    N(i).spktrain = N(i).spktrain(1:spT,:);
    S(i).spktrain = S(i).spktrain(1:spT,:);
    
    N(i).spktrain_bl = N(i).spktrain_bl(1:spT,:);
    S(i).spktrain_bl = S(i).spktrain_bl(1:spT,:);
end

Nspk = [cell2mat({N.spktrain_bl}); cell2mat({N.spktrain})];
Sspk = [cell2mat({S.spktrain_bl}); cell2mat({S.spktrain})];
Fs = 10000;
bin = 111.1;
[spbinN] = psth_sp(Nspk,bin,Fs);
[spbinS] = psth_sp(Sspk,bin,Fs);

spbinN = smooth(spbinN);
spbinS = smooth(spbinS);
maxfr = max(spbinN);
spbinN = spbinN ./ maxfr;
spbinS = spbinS ./ maxfr;
t = linspace(-spT/10,spT/10,size(spbinN,1));

NN = spbinN;
SS = spbinS;
NN(51) = NaN;
SS(51) = NaN;

h1 = figure;
% subplot(2,1,2) 
plot(t,NN,'k','LineWidth',2); hold on
plot(t,SS,'Color',[0 .5 0],'LineWidth',2) % t(1:end-1)
% line([150 150],[0 30],'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2,'LineStyle','--')
shade(t(50:52),ones(3,1)*ceil(max(NN)),'FillColor',[.5 .5 .5])
xlim([-400 400])
ylim([0 ceil(max(NN))])
box off
title('MT Population','FontSize',15,'FontName','High Tower Text')
xlabel('Time(sec)','FontSize',15,'FontName','High Tower Text')
ylabel('Firing Rate_{p.u.}','FontSize',15,'FontName','High Tower Text')
legend({'Control','MS'},'box','off')
h1.Color = 'w';
[hstim,pstim] = ttest(NN(37:50),SS(37:50))
[hstimE,pstimE] = ttest(NN(52:65),SS(52:65))
stim_avg = mean(SS(37:50)) - mean(NN(37:50))
stimE_avg = mean(SS(52:65)) - mean(NN(52:65))

%%
for i = 13
celln = i;
an = horzcat(Ntrans.info.sp);
bn =  horzcat(Nspiral.info.sp);
cn = [an(celln) bn(celln)];
cellN = vertcat(cn(1).spktrain_bl(end-4050:end,:),cn(1).spktrain);
cellN = [cellN vertcat(cn(2).spktrain_bl(end-4050:end,:),cn(2).spktrain)];

as = horzcat(Strans.info.sp);
bs =  horzcat(Sspiral.info.sp);
cs = [as(celln) bs(celln)];
cellS = vertcat(cs(1).spktrain_bl(end-4050:end,:),cs(1).spktrain);
cellS = [cellS vertcat(cs(2).spktrain_bl(end-4050:end,:),cs(2).spktrain)];

pstn = psth_sp(cellN,bin,Fs);
psts = psth_sp(cellS,bin,Fs);

pstn = smooth(pstn);
psts = smooth(psts);
pstn(51) = NaN;
psts(51) = NaN;

h2 = figure;
% subplot(2,1,1)
plot(t,pstn,'k','LineWidth',2); hold on
plot(t,psts,'Color',[0 .5 0],'LineWidth',2)
shade(t(50:52),ones(3,1)*ceil(max(pstn)),'FillColor',[.5 .5 .5])
xlim([-400 400])
h2.Color = 'w';
box off
title('Neuron #13','FontSize',15,'FontName','High Tower Text')
% set(gca,'xTickLabel','')
% ylim([0 ceil(max(pstn))])
xlabel('Time(sec)','FontSize',15,'FontName','High Tower Text')
ylabel('Firing Rate','FontSize',15,'FontName','High Tower Text')
legend({'Control','MS'},'box','off')


end
%%
Nspk = [cell2mat({N.spktrain_bl}); cell2mat({N.spktrain})];
Sspk = [cell2mat({S.spktrain_bl}); cell2mat({S.spktrain})];
Fs = 10000;
bin = 1;
[spbinN] = psth_sp(Nspk,bin,Fs);
[spbinS] = psth_sp(Sspk,bin,Fs);
t = linspace(-spT/10,spT/10,size(spbinN,1));
NN = spbinN;
SS = spbinS;
% NN(5513:5738) = NaN;
% SS(5513:5738) = NaN;
% figure
subplot(2,1,1)
plot(t,NN,'k','LineWidth',2); hold on
plot(t,SS,'Color',[0 .5 0],'LineWidth',2)
uylimit = max(spbinS);
% line([150 150],[0 uylimit],'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2,'LineStyle','--')
shade(t(5512:5739),ones(length(t(5512:5739)),1)*600,'FillColor',[.5 .5 .5])
xlim([-400 400])
ylim([0 600])
set(gca,'xTickLabel','')
h.Color = 'w';
box off
ylabel('Spikes/sec','FontSize',15,'FontName','High Tower Text')

bar(t,NN,50,'BarWidth',1,'FaceColor','k','EdgeColor','k'); hold on
bar(t,SS,50,'BarWidth',1,'FaceColor',[0 .5 0],'EdgeColor',[0 .5 0])
alpha(.1)
legend({'Control','MS'},'box','off')
