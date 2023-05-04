clear all; clc; %close all
%% This is just to make a grid and then superpose an orientation map on the grid 
grid = 81;    % grid size
netpos = linspace(-5,5,grid);
%% Compute the distance between grid intervals: regardless of E/I, the distance difference of each to grid interval is stored here
deltax = deg2rad(10/(grid-1));
[X,Y] = meshgrid(reshape(netpos',1,grid),reshape(netpos',1,grid));
xyD = distmat(X');
xyD = xyD * deltax; %scale to deg
% figure; imagesc(xyD); axis square
% set(gca,'YDir','Normal')
%% parameters
Jee = 0.385;
Jie = 1;
Jei = 0.55;
Jii = 0.018;
sigmaRF = .33 * deltax;
kappaE = .1;
kappaIE = 1;
SigmaEE = deg2rad(0.5);
SigmaIE = deg2rad(1); 
%% Compute Weights
GsigmaDirEE = kappaE * exp(-(xyD.^2) / (2 * SigmaEE^2));
GsigmaDirIE = kappaIE * exp(-(xyD.^2) / (2 * SigmaIE^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wee = rand(grid,grid);
Wee(Wee < GsigmaDirEE) = 1;
Wee(Wee ~= 1) = 0;
Wee = Wee .* normrnd(Jee,0.25*Jee,[grid grid]);
Wee = max(Wee,0);
cWee = Jee * sum(GsigmaDirEE,2);
aWee = sum(Wee,2);
sWee = cWee ./ aWee;
sWee(isinf(sWee)==1) = 0;
Wee = Wee .* repmat(sWee,1,grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Wie = rand(grid,grid);
Wie(Wie < GsigmaDirIE) = 1;
Wie(Wie ~= 1) = 0;
Wie = Wie .* normrnd(Jie,0.25*Jie,[grid grid]);
Wie = max(Wie,0);
cWie = Jie * sum(GsigmaDirIE,2);
aWie = sum(Wie,2);
sWie = cWie ./ aWie;
sWie(isinf(sWie)==1) = 0;
Wie = Wie .* repmat(sWie,1,grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Wii = zeros(grid,grid);
Wei = zeros(grid,grid);
Wii(find(eye(size(Wii)))) = 1.5;
Wei(find(eye(size(Wei)))) = 0.55;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
W = [Wee -Wei; Wie -Wii];
figure
subplot(2,2,1); imagesc(Wee); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
title([sprintf('%.4f',min(Wee(:))),'\leq','W_{EE}','\leq',sprintf('%.4f',max(Wee(:)))])
subplot(2,2,2); imagesc(Wei); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
title([sprintf('%.4f',min(Wei(:))),'\leq','W_{EI}','\leq',sprintf('%.4f',max(Wei(:)))])
subplot(2,2,3); imagesc(Wie); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
title([sprintf('%.4f',min(Wie(:))),'\leq','W_{IE}','\leq',sprintf('%.4f',max(Wie(:)))])
subplot(2,2,4); imagesc(Wii); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
title([sprintf('%.4f',min(Wii(:))),'\leq','W_{II}','\leq',sprintf('%.4f',max(Wii(:)))])
colormap(copper)
xlabel '-5^o,...,0^o,..,5^o'
ylabel '-5^o,...,0^o,..,5^o'
suptitle 'Connection Weights'
%% h(theta)
cen = round(grid/2);
lrange = linspace(0.1,10,20);
l = deg2rad(lrange);
X = 1:grid;
X = X - cen;
X = X * deltax;
% figure
hpos = [];
for i = 1:size(lrange,2)
    hpos(i,:) = (1 ./ (1 + exp(-((X + (l(i)/2)) / sigmaRF))))...
         + (1 - (1 ./ (1 + exp(-((X - (l(i)/2)) / sigmaRF))))) - 1;
%     plot(netpos,hpos(i,:)); hold on
end

X = 1:grid;
X = X - cen-23;
X = X * deltax;
hp = (1 ./ (1 + exp(-((X + (l(6)/2)) / sigmaRF))))...
         + (1 - (1 ./ (1 + exp(-((X - (l(6)/2)) / sigmaRF))))) - 1;
%% Tau, k, & n
makeitRND = 'no';
TSpan = [0 400];
r0 = zeros(2*grid,1);
tauE = 20;
tauI = 10;
tau = zeros(2*grid,1);
switch makeitRND
    case 'yes'
        tau = [normrnd(tauE,0.05 * tauE, [grid 1]); normrnd(tauI,0.05 * tauI, [grid 1])];
        unE = 2.2;
        unI = 2.2;
        n = [normrnd(unE,0.05*unE, [grid 1]);normrnd(unI,0.05*unI, [grid 1])];
        uk = 0.01;
        k = normrnd(uk,0.05*uk, [2*grid 1]);
    case 'no'
        tau = [normrnd(tauE,0.05 * tauE, [grid 1]); normrnd(tauI,0.05 * tauI, [grid 1])];
        n = ones(2*grid,1);
        k = ones(2*grid,1);
end

%%
C = [1;0;1;3;1];
A = [0;1;1;1;2];
for j = 1:5
    RE = []; RI = []; Fr = [];
    c = C(j);
    a = A(j);
    for i = 1:20
        We = [Wee -Wei; zeros(grid,grid) zeros(grid,grid)];
        [te,re] = ode45(@(te,re) isn2E(te,re,sparse(We),[hpos(i,:)';zeros(81,1)],tau,c,k,n,a,[hp';zeros(81,1)]),TSpan,r0);
        Fre = real(re(end,:));
        RE(:,i) = Fre(1:grid);
        
        Wi = [zeros(grid,grid) zeros(grid,grid); Wie -Wii];
        [ti,ri] = ode45(@(ti,ri) isn2(ti,ri,sparse(Wi),[zeros(81,1);hpos(i,:)'],tau,1,k,n),TSpan,r0);
        Fri = real(ri(end,:));
        RI(:,i) = Fri(grid+1:end);
        
        h = [hpos(i,:)';hpos(i,:)'];
        [t,r] = ode45(@(t,r) isn2(t,r,sparse(W),h,tau,1,k,n),TSpan,r0);
        Fr(:,i) = real(r(end,:));
    end
    x{j,1} = RE;
    x{j,2} = RI;
    x{j,3} = Fr;
end

%%
stim = 6;
col = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0, 0.5, 0];[0.6350, 0.0780, 0.1840];[0.4940, 0.1840, 0.5560]];
lg = {'C';'S';'C+S';'C+S,Attend C';'C+S,Attend S'};
figure;
subplot(6,1,1)
plot(netpos,hpos(stim,:)','Color',col(1,:),'LineWidth',1.5); hold on
plot(netpos,hp','Color',col(2,:),'LineWidth',1.5)
set(gca,'Box','off')
for i = 1:5
    subplot(6,1,i+1)
    plot(netpos,x{i,1}(:,stim)','Color',col(i,:),'LineWidth',1.5)
    set(gca,'Box','off')
%     axis square
    legend(lg{i},'Location','northoutside')
end
xlabel 'Network Position'
ylabel 'Firing Rate'
set(gcf,'color','w')
%%
rE = RE;
rI = RI;
figure;
for i = 1:size(lrange,2)/2
    subplot(3,size(lrange,2)/2,i)
    plot(netpos,hpos(i,:),'k','LineWidth',1)
    set(gca,'Box','off')
    ylim([0 max(hpos(:))+0.05])
    if i == 1
        ylabel 'Input'
    end
    
    subplot(3,size(lrange,2)/2,i+size(lrange,2)/2)
    plot(netpos,rE(:,i)','r','LineWidth',.5)
    set(gca,'Box','off')
    ylim([min(rE(:))-0.01 max(rE(:))+0.01])
    
    subplot(3,size(lrange,2)/2,i+(size(lrange,2)/2*2))
    plot(netpos,rI(:,i)','b','LineWidth',.5)
    set(gca,'Box','off')
    ylim([min(rI(:))-0.01 max(rI(:))+0.01])
end
axE = subplot(3,size(lrange,2)/2,11);
axE.YLabel.String = 'Firing Rate';
axI = subplot(3,size(lrange,2)/2,21);
axI.YLabel.String = axE.YLabel.String;
axI.XLabel.String = 'Network Position';
%%
% figure
% subplot(3,1,1)
% plot(hpos(i,:),'k')
% 
% subplot(3,1,2)
% plot(Fre(1:81),'r')
% hold on
% plot(Fri(82:end),'b')
% 
% subplot(3,1,3)
% plot(Fr(:,1:grid),'r')
% hold on
% plot(Fr(:,grid+1:grid*2),'b')
%% Solve ODE
% c = 1;
% for i = 1:size(lrange,2)
%     h = [hpos(i,:)';hpos(i,:)'];
%     [t,r] = ode45(@(t,r) isn2(t,r,sparse(W),h,tau,c,k,n),TSpan,r0);
%     Fr = real(r(end,:));
%     rE(:,i) = Fr(:,1:grid);
%     rI(:,i) = Fr(:,grid+1:grid*2);
%     rEc(1,i) = Fr(ceil(grid/2));
%     rIc(1,i) = Fr(ceil(grid/2)*3);
% end
%% Figures
rE = RE;%Fr(1:grid,:);
rI = RI;%Fr(grid+1:end,:);
figure;
for i = 1:size(lrange,2)/2
    subplot(3,size(lrange,2)/2,i)
    plot(netpos,hpos(i,:),'k','LineWidth',1)
    set(gca,'Box','off')
    ylim([0 max(hpos(:))+0.05])
    if i == 1
        ylabel 'Input'
    end
    
    subplot(3,size(lrange,2)/2,i+size(lrange,2)/2)
    plot(netpos,rE(:,i)','r','LineWidth',.5)
    set(gca,'Box','off')
    ylim([min(rE(:))-0.01 max(rE(:))+0.01])
    
    subplot(3,size(lrange,2)/2,i+(size(lrange,2)/2*2))
    plot(netpos,rI(:,i)','b','LineWidth',.5)
    set(gca,'Box','off')
    ylim([min(rI(:))-0.01 max(rI(:))+0.01])
end
axE = subplot(3,size(lrange,2)/2,11);
axE.YLabel.String = 'Firing Rate';
axI = subplot(3,size(lrange,2)/2,21);
axI.YLabel.String = axE.YLabel.String;
axI.XLabel.String = 'Network Position';
return
%%
figure;
subplot(2,1,1)
plot(rEc,'--r*','LineWidth',1)
ylabel 'Firing Rate'
xlim([1 20])
set(gca,'Box','off')
subplot(2,1,2)
plot(rIc,'--b*','LineWidth',1)
xlim([1 20])
xlabel 'Stimulus Length'
ylabel 'Firing Rate'
set(gca,'Box','off')