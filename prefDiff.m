function [pdDiff Fr] = prefDiff(datapath,bin)
%%
dp = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
Nt = load([dp 'Ntrans.mat']);
Ns = load([dp 'Nspiral.mat']);
St = load([dp 'Strans.mat']);
Ss = load([dp 'Sspiral.mat']);
mstNt = load([dp 'NtransMS.mat']);
mstNs = load([dp 'NspiralMS.mat']);


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

tr = struct2cell(mstNt);
sp = struct2cell(mstNs);
for i = 1:size(tr{1,1},2)
    for j = 1:size(tr{1,1}(i).xr,2)
        tr{1,1}(i).xr(j).firing = tr{1,1}(i).xr(j).firing./max(tr{1,1}(i).xr(j).firing(:));
        sp{1,1}(i).xr(j).firing = sp{1,1}(i).xr(j).firing./max(sp{1,1}(i).xr(j).firing(:));
    end
end
NT = [tr{1,1}(:).xr];
NS = [sp{1,1}(:).xr];

%%%% Input Data Section %%%%
load([datapath 'theta_t_MS.mat']); trMSTt = theta; clear theta
load([datapath 'theta_s_MS.mat']); trMSTs = theta; clear theta
% load([datapath 'trMST_t_MS.mat']); trMSTt = trMST; clear theta
% load([datapath 'trMST_s_MS.mat']); trMSTs = trMST; clear theta
load([datapath 'theta_t_n.mat']); trMTt = theta; clear theta; frMTt = mfr; clear mfr
load([datapath 'theta_s_n.mat']); trMTs = theta; clear theta; frMTs = mfr; clear mfr

%%%% Computation Section %%%%
c = 1;
for i = 1:max(trMTt(:,7))
    mtidx = find(trMTt(:,7)==i);
    %     mstT = trMSTt(i);
    %     mstS = trMSTs(i);
    mstT = trMSTt(trMSTt(:,4)==i,1);
    mstTFr = trMSTt(trMSTt(:,4)==i,3);
    mstS = trMSTs(trMSTs(:,4)==i,1);
    mstSFr = trMSTs(trMSTs(:,4)==i,3);
    mstRFt = tr{1,1}(i).xr;
    mstRFs = sp{1,1}(i).xr;
    mtRFt = N1{1,1}(i).xr;
    mtRFs = N2{1,1}(i).xr;
    for j = 1:size(mstT,1)
        for z = 1:size(mtidx,1)
            pdDiff(c,1) = abs(angdiff(trMTt(mtidx(z),1),mstT(j))); % MST Trans NoMS vs. MT Trans MS
            pdDiff(c,2) = abs(angdiff(trMTt(mtidx(z),2),mstT(j))); % MST Trans NoMS vs. MT Trans NoMS
            pdDiff(c,3) = abs(angdiff(trMTs(mtidx(z),1),mstS(j))); % MST Spiral NoMS vs. MT Spiral MS
            pdDiff(c,4) = abs(angdiff(trMTs(mtidx(z),2),mstS(j))); % MST Spiral NoMS vs. MT Spiral NoMS
            [Vxt(c,:,:) Vyt(c,:,:) diffVt(c,:,:)] = RFupdated(mstRFt(j),mtRFt(z));
            [Vxs(c,:,:) Vys(c,:,:) diffVs(c,:,:)] = RFupdated(mstRFs(j),mtRFs(z));
            RF(c,1) = abs(Edist(Vxt(c,:,:),Vyt(c,:,:)));
            RF(c,2) = abs(Edist(Vxs(c,:,:),Vys(c,:,:)));
            %             Fr(c,1) = mstTFr(j,1) - trMTt(mtidx(z),5); % MS
            %             Fr(c,2) = mstTFr(j,1) - trMTt(mtidx(z),6); % NMS
            %             Fr(c,3) = mstSFr(j,1) - trMTs(mtidx(z),5); % MSclose all
            %             Fr(c,4) = mstSFr(j,1) - trMTs(mtidx(z),6); % NMS
            Fr(c,1) = frMTt(mtidx(z),1); % MS
            Fr(c,2) = frMTt(mtidx(z),2); % NMS
            Fr(c,3) = frMTs(mtidx(z),1); % MS
            Fr(c,4) = frMTs(mtidx(z),2); % NMS
            %             Fr(c,1) = trMTt(mtidx(z),5); % MS
            %             Fr(c,2) = trMTt(mtidx(z),6); % NMS
            %             Fr(c,3) = trMTs(mtidx(z),5); % MS
            %             Fr(c,4) = trMTs(mtidx(z),6); % NMS
            Fr(c,5) = mstTFr(j,1); % MST NMS Trans
            Fr(c,6) = mstSFr(j,1); % MST NMS Spiral
            Fr(c,7) = i;
            c = c + 1;
        end
    end
end
%% Delta FR vs Delta PD
frt = Fr(:,1) - Fr(:,2);
frs = Fr(:,3) - Fr(:,4);
supermin = min(min(frs),min(frt));
supermax = max(max(frs),max(frt));
offset = 0.1;

x = [];y=[];X=[];Y=[];H=[];
x = vertcat(pdDiff(:,2),pdDiff(:,4));
y = vertcat(frt,frs);
X = [ones(length(x),1) x];
H = X\y;
Y = X * H;
mdl1 = fitlm(x,y);
Rsq1 = mdl1.Rsquared.Ordinary;
Pval1 = mdl1.Coefficients.pValue(2);
[pSign1,hSign1,stat1] = signrank(y);
% Plot dFR vs dPD
h1 = figure;
patch([0 pi+0.1 pi+0.1 0],[0 0 supermin-offset supermin-offset],[0.5 0.5 0.5],'FaceAlpha',.4,'EdgeColor',[.94 .94 .94])
patch([0 pi+0.1 pi+0.1 0],[0 0 supermax+offset supermax+offset],[0.5 0.5 0.5],'FaceAlpha',.2,'EdgeColor',[.94 .94 .94])
hold on
scatter(pdDiff(:,2),frt,10,'k','filled')
hold on
scatter(pdDiff(:,4),frs,10,'k','filled')
plot(x,Y,'LineWidth',2,'Color',[0, 0.5, 0])
h1.Color = 'w';
xlim([0 pi+0.1])
ylim([supermin-offset supermax+offset])
% Modify labels based on MST-MT or MST-MST input
xlabel('\DeltaPd_{Control}^{r}','FontSize',15,'FontName','High Tower Text')
ylabel('\DeltaFr_{p.u.}','FontSize',15,'FontName','High Tower Text')
title(sprintf('R^2=%.2f p-val=%.3f',Rsq1,Pval1))
h1_1 = figure;
histogram(y,30,'FaceColor','[0.75, 0.75, 0]','EdgeColor','[0.75, 0.75, 0]','Normalization','count','FaceAlpha',.8)
xlim([min(y)-0.01 max(y)+0.01])
box off
h1_1.Color = 'w';
h1_1.CurrentAxes.YAxis.Visible = 'off';
% Plot dFR vs dPD by motion type
% figure
% subplot(1,2,1)
% patch([0 pi+0.1 pi+0.1 0],[0 0 supermin-offset supermin-offset],'b','FaceAlpha',.2)
% patch([0 pi+0.1 pi+0.1 0],[0 0 supermax+offset supermax+offset],'r','FaceAlpha',.2)
% hold on
% scatter(pdDiff(:,2),frt,'b','filled')
% xlim([0 pi+0.1])
% ylim([supermin-offset supermax+offset])
% title 'Translation'
% xlabel '{\Delta}PD^o = {\Theta}_{MST} - {\Theta}_{MST}^{NMS}'
% ylabel '{\Delta}FR (spk/sec) = Fr_{MT}^{MS} - Fr_{MST}^{NMS}'
% 
% subplot(1,2,2)
% patch([0 pi+0.1 pi+0.1 0],[0 0 supermin-offset supermin-offset],'b','FaceAlpha',.2)
% patch([0 pi+0.1 pi+0.1 0],[0 0 supermax+offset supermax+offset],'r','FaceAlpha',.2)
% hold on
% scatter(pdDiff(:,4),frs,'k','filled')
% xlim([0 pi+0.1])
% ylim([supermin-offset supermax+offset])
% title 'Spiral'
% xlabel '{\Delta}PD^o = {\Theta}_{MST} - {\Theta}_{MST}^{NMS}'
%% Delta PDms vs Delta PDnms
x = [];y=[];X=[];Y=[];H=[];
x = vertcat(pdDiff(:,2),pdDiff(:,4));
y = vertcat(pdDiff(:,1),pdDiff(:,3));
X = [ones(length(x),1) x];
H = X\y;
Y = X * H;
mdl2 = fitlm(x,y);
Rsq2 = mdl2.Rsquared.Ordinary;
Pval2 = mdl2.Coefficients.pValue(2);
% [pSign2,hSign2,stat2] = signrank(x,y);

h2 = figure;
oset = 3.2;
patch([0 oset oset 0],[0 0 oset 0],[0.5 0.5 0.5],'FaceAlpha',.4,'EdgeColor',[.94 .94 .94])
patch([0 oset 0 0],[0 oset oset 0],[0.5 0.5 0.5],'FaceAlpha',.2,'EdgeColor',[.94 .94 .94])
hold on
scatter(pdDiff(:,2),pdDiff(:,1),10,'MarkerEdgeColor','k','MarkerFaceColor','k')
scatter(pdDiff(:,4),pdDiff(:,3),10,'MarkerEdgeColor','k','MarkerFaceColor','k')
plot(x,Y,'LineWidth',2,'Color',[0, 0.5, 0])
box off
xlim([0 oset]); ylim([0 oset])
% Modify labels based on MST-MT or MST-MST input
xlabel('\DeltaPd_{Control}^{r}','FontSize',15,'FontName','High Tower Text')
ylabel('\DeltaPd_{Microstim}^{r}','FontSize',15,'FontName','High Tower Text')
h2.Color = 'w';
title(sprintf('R^2=%.2f p-val=%.3f',Rsq2,Pval2))
%%%% Plot Angle Difference Distribution %%%%
h2_1 = figure;
z = angdiff(y,x);
[pSign2_1,hSign2_1,stat2_1] = signrank(z)
histogram(z,30,'FaceColor','[0.75, 0.75, 0]','EdgeColor','[0.75, 0.75, 0]','Normalization','count','FaceAlpha',.8)
box off
% xlabel('Difference^{rad}','FontSize',15,'FontName','High Tower Text')
% ylabel('Counts','FontSize',15,'FontName','High Tower Text')
h2_1.Color = 'w';
xlim([-oset oset])
h2_1.CurrentAxes.YAxis.Visible = 'off'; 
%% Delta FR vs Delta RF
fr = [frt;frs];
rf = [RF(:,1);RF(:,2)];

x = [];y=[];X=[];Y=[];H=[];
x = rf;
y = fr;
X = [ones(length(x),1) x];
H = X\y;
Y = X * H;
mdl3 = fitlm(x,y);
Rsq3 = mdl3.Rsquared.Ordinary;
Pval3 = mdl3.Coefficients.pValue(2);
[pSign3,hSign3,stat3] = signrank(y);

h3 = figure;
patch([0 max(rf)+0.1 max(rf)+0.1 0],[0 0 supermin-offset supermin-offset],[0.5 0.5 0.5],'FaceAlpha',.4,'EdgeColor',[.94 .94 .94])
patch([0 max(rf)+0.1 max(rf)+0.1 0],[0 0 supermax+offset supermax+offset],[0.5 0.5 0.5],'FaceAlpha',.2,'EdgeColor',[.94 .94 .94])
hold on
scatter(rf,fr,10,'k','filled')
plot(x,Y,'LineWidth',2,'Color',[0, 0.5, 0])
h3.Color = 'w';
xlim([0 max(rf)+0.1])
ylim([supermin-offset supermax+offset])
% title 'MT Firing Rate across {\Delta}PD'
xlabel('|\DeltaRF_{Control}|','FontSize',15,'FontName','High Tower Text')
ylabel('\DeltaFr_{p.u.}','FontSize',15,'FontName','High Tower Text')
title(sprintf('R^2=%.2f p-val=%.3f',Rsq3,Pval3))
%%
fr = [];
fr = [frMTt;frMTs];
pwpNMS = pwproduct(fr(:,2));
pwpMS = pwproduct(fr(:,1));
% figure
% scatter(pwpNMS,pwpMS)
figure
pwp = pwpMS - pwpNMS;
hist(pwp,50)
ttest(pwp)
%% For RF and other extras
% dp = 'D:\MT_MST\Microstim\Norm\MUA_SUA\';
% Nt = load([dp 'Ntrans.mat']);
% Ns = load([dp 'Nspiral.mat']);
% St = load([dp 'Strans.mat']);
% Ss = load([dp 'Sspiral.mat']);
% mstNt = load([dp 'NtransMS.mat']);
% mstNs = load([dp 'NspiralMS.mat']);
% 
% 
% Nt.info(7) = []; Nt.info(4) = []; Nt.info(1) = [];
% Ns.info(7) = []; Ns.info(4) = []; Ns.info(1) = [];
% 
% St.info(7) = []; St.info(4) = []; St.info(1) = [];
% Ss.info(7) = []; Ss.info(4) = []; Ss.info(1) = [];
% 
% N1 = struct2cell(Nt);
% N2 = struct2cell(Ns);
% Nt = [N1{1,1}(:).xr];
% Ns = [N2{1,1}(:).xr];
% 
% S1 = struct2cell(St);
% S2 = struct2cell(Ss);
% St = [S1{1,1}(:).xr];
% Ss = [S2{1,1}(:).xr];
% 
% tr = struct2cell(mstNt);
% sp = struct2cell(mstNs);
% NT = [tr{1,1}(:).xr];
% NS = [sp{1,1}(:).xr];
% 
% [Vxt Vyt diffVt] = RFupdated(St,Nt);
% [Vxs Vys diffVs] = RFupdated(Ss,Ns);
% 
% for i = 1:size(Vxt,1)
%     RF(i,1) = abs(Edist(Vxt(i,:),Vyt(i,:))); % Translation
%     RF(i,2) = abs(Edist(Vxs(i,:),Vys(i,:))); % Spirals
% end
% %%
% FrRf = zeros(size(frMTt,1)*2,2);
% FrRf(:,1) = [(frMTt(:,1) - frMTt(:,2));(frMTs(:,1) - frMTs(:,2))];
% FrRf(:,2) = [RF(:,1);RF(:,2)];
% FrRf(44,:) = [];
% %%
% x = [];y=[];X=[];Y=[];H=[];
% x = FrRf(:,2);
% y = FrRf(:,1);
% X = [ones(length(x),1) x];
% H = X\y;
% Y = X * H;
% 
% rff = figure;
% scatter(FrRf(:,2),FrRf(:,1),10,'MarkerEdgeColor','k','MarkerFaceColor','k')
% hold on
% plot(x,Y,'LineWidth',2,'Color',[.5 .5 .5])
% box off
% xlabel('|\DeltaRF|_{MS-NoMs}^{\circ}','FontSize',15,'FontName','High Tower Text')
% ylabel('\DeltaFR_{MS-NoMs}','FontSize',15,'FontName','High Tower Text')
% rff.Color = 'w';

%%%% Plot CDF %%%%
% figure
% subplot(2,2,1)
% F1 = cdfplot(pdDiff(:,1));
% set(F1,'LineWidth',2,'Color','b')
% hold on
% F2 = cdfplot(pdDiff(:,2));
% set(F2,'LineWidth',2,'Color','g')
% legend({'MS','NMS'},'Location','East')
% grid off
% title 'CDF PD Diff'
% xlabel ''
% ylabel 'Probability'
% 
% subplot(2,2,2)
% F1 = cdfplot(pdDiff(:,3));
% set(F1,'LineWidth',2,'Color','k')
% hold on
% F2 = cdfplot(pdDiff(:,4));
% set(F2,'LineWidth',2,'Color','r')
% legend({'MS','NMS'},'Location','East')
% grid off
% title 'CDF PD Diff'
% xlabel ''
% ylabel ''
% 
% subplot(2,2,3)
% histogram(pdDiff(:,1),bin,'FaceColor','b')
% hold on
% histogram(pdDiff(:,2),bin,'FaceColor','g')
% title(sprintf('PD_{Diff}^{Trans} p-val = %.2f',pt))
% xlabel('(MST - MT)^o_{Diff}')
% ylabel('#of pairs')
% legend({'MS','NMS'})
% xlim([0 3.15])
% 
% subplot(2,2,4)
% histogram(pdDiff(:,3),bin,'FaceColor','k')
% hold on
% histogram(pdDiff(:,4),bin,'FaceColor','r')
% title(sprintf('PD_{Diff}^{Spiral} p-val = %.2f',ps))
% legend({'MS','NMS'})
% xlim([0 3.15])
%%%% Plot dPDms vs dPDnms based on motion type %%%%
% figure
% subplot(1,2,1)
% scatter(pdDiff(:,2),pdDiff(:,1),'MarkerEdgeColor','g','MarkerFaceColor','b'); refline(1,0)
% xlim([0 3.15]); ylim([0 3.15])
% title('PD_{Diff}^{Trans}')
% xlabel('(MST_{NoMS} - MT_{NoMS})^o')
% ylabel('(MST_{NoMS} - MT_{MS})^o')
% 
% subplot(1,2,2)
% scatter(pdDiff(:,4),pdDiff(:,3),'MarkerEdgeColor','r','MarkerFaceColor','k'); refline(1,0)
% xlim([0 3.15]); ylim([0 3.15])
% title('PD_{Diff}^{Spiral}')
% xlabel('(MST_{NoMS} - MT_{NoMS})^o')
% pdD = zeros(size(pdDiff,1),2);
% pdD(:,1) = abs(angdiff(pdDiff(:,1),pdDiff(:,2)));
% pdD(:,2) = abs(angdiff(pdDiff(:,3),pdDiff(:,4)));
% [h,p] = kstest2(pdD(:,1),pdD(:,2));
% 
% figure
% subplot(2,1,1)
% F1 = cdfplot(pdD(:,1));
% set(F1,'LineWidth',2,'Color','b')
% hold on
% F2 = cdfplot(pdD(:,2));
% set(F2,'LineWidth',2,'Color','k')
% legend('Trans','Spiral','Location','East')
% grid off
% title 'CDF PD'
% xlabel ''
% ylabel 'Probability'
% 
% subplot(2,1,2)
% histogram(pdD(:,1),bin,'FaceColor','b')
% hold on
% histogram(pdD(:,2),bin,'FaceColor','k')
% xlim([0 pi])
% % legend('Trans','Spiral')
% title(sprintf('MT^{PD} Diff^{NoMS - MS} from MST^{PD} p-val = %.2f',pt))
% xlabel('Degree (rad)')
% ylabel('Paris of MT-MST')
