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
% figure
% subplot(2,2,1); imagesc(Wee); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
% title([sprintf('%.4f',min(Wee(:))),'\leq','W_{EE}','\leq',sprintf('%.4f',max(Wee(:)))])
% subplot(2,2,2); imagesc(Wei); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
% title([sprintf('%.4f',min(Wei(:))),'\leq','W_{EI}','\leq',sprintf('%.4f',max(Wei(:)))])
% subplot(2,2,3); imagesc(Wie); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
% title([sprintf('%.4f',min(Wie(:))),'\leq','W_{IE}','\leq',sprintf('%.4f',max(Wie(:)))])
% subplot(2,2,4); imagesc(Wii); axis square; set(gca,'YDir','Normal','YTick',[],'XTick',[])
% title([sprintf('%.4f',min(Wii(:))),'\leq','W_{II}','\leq',sprintf('%.4f',max(Wii(:)))])
% colormap(copper)
% xlabel '-5^o,...,0^o,..,5^o'
% ylabel '-5^o,...,0^o,..,5^o'
% suptitle 'Connection Weights'
%% h(theta)
cen = round(grid/2);
lrange = linspace(0.1,10,20);
l = deg2rad(lrange);
X = 1:grid;
X = X - cen;
X = X * deltax;
stimsize = 6;
hpos = (1 ./ (1 + exp(-((X + (l(stimsize)/2)) / sigmaRF))))...
         + (1 - (1 ./ (1 + exp(-((X - (l(stimsize)/2)) / sigmaRF))))) - 1;
% FB
x = 1:grid;
m1 = cen - floor(cen/2);
m2 = cen + floor(cen/2);
sigma1 = 0.01;
sigma2 = 0.01;
% FBe = bimodal_dist(x,m1,m2,sigma1,sigma2,.01,.01);
% FBi = bimodal_dist(x,m1,m2,sigma1,sigma2,.02,.02);

e = 10;
i = 12;

FBe = (1 ./ (1 + exp(-((X + (l(e)/2)) / sigma1))))...
    + (1 - (1 ./ (1 + exp(-((X - (l(e)/2)) / sigma1))))) - 1;

FBi = (1 ./ (1 + exp(-((X + (l(i)/2)) / sigma2))))...
    + (1 - (1 ./ (1 + exp(-((X - (l(i)/2)) / sigma2))))) - 1;

figure
plot(hpos,'k','LineWidth',2); hold on
plot(FBe,'r','LineWidth',2)
plot(FBi,'b','LineWidth',2)
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
        unE = 2.0;
        unI = 2.2;
        n = [normrnd(unE,0.05*unE, [grid 1]);normrnd(unI,0.05*unI, [grid 1])];
        uk = 0.01;
        k = normrnd(uk,0.05*uk, [2*grid 1]);
    case 'no'
        tau = [normrnd(tauE,0.05 * tauE, [grid 1]); normrnd(tauI,0.05 * tauI, [grid 1])];
        n = ones(2*grid,1);
        k = ones(2*grid,1);
end

%% Solve ODE
FB = [-FBe,-FBi;zeros(1,grid),-FBi;-FBe,zeros(1,grid)];
a = [0;1];

figure
b = [{'FB_E + FB_I'},{'0_E + FB_I'},{'FB_E + 0_I'}];
subplot(2,1,1)
plot(netpos,hpos,'k','LineWidth',2); hold on
plot(netpos,FBe,'r','LineWidth',2)
plot(netpos,FBi,'b','LineWidth',2)
legend({'h(x)','E','I'})
xlim([-5 5])
xlabel 'x_{-5^o,...,0^o,..,5^o}'
set(gca,'Box','off')

for j = 1%:3
    for i = 1:2
        [t,r] = ode45(@(t,r) isn2E(t,r,sparse(W),[hpos';hpos'],tau,1,k,n,a(i),FB(j,:)'),TSpan,r0);
        Fr(j,i,:) = real(r(end,:));
        RE(j,i,:) = Fr(j,i,1:grid);
        RI(j,i,:) = Fr(j,i,grid+1:end);
        R = squeeze([RE(j,i,:);RE(j,i,:)]);
        R = mean(R);
        R1 = flip(R(1:cen));
        R2 = R(cen:end);
        Rmean(j,i,:) = mean([R1;R2]);
    end
    subplot(2,1,j+1)
    plot(netpos(cen:end),squeeze(Rmean(j,1,:)),'k','LineWidth',1.5); hold on
    plot(netpos(cen:end),squeeze(Rmean(j,2,:)),'b','LineWidth',1.5)
    title(b{j})
    set(gca,'Box','off')
end
xlabel 'Dev. from Network Center'
ylabel 'Mean Firing Rate of E/I pair'
legend({'NMS','MS'})
set(gcf,'color','w')


