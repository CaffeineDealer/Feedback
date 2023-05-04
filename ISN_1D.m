clear all; clc; %close all
%% This is just to make a grid and then superpose an orientation map on the grid 
grid = 360;    % grid size
th_size = 360;
theta = zeros(1,th_size);

theta = deg2rad(linspace(0,359,th_size));
figure; 
imagesc(theta); axis square
colormap(bone); colorbar
title 'Map of PD for each E/I pair'
xlabel 'x_{0^o,1^o,..,359^o}'
set(gca,'YTick',[],'XTick',[])
%% parameters
deltatheta = deg2rad(1);
Jee = 0.044;
Jie = 0.042;
Jei = 0.023;
Jii = 0.018;
sigmadir = deg2rad(64);
sigmaFF = deg2rad(60);
kappaE = 0.1;
kappaI = 0.2;
%% Compute Direction preference for each pair
[X,Y] = meshgrid(reshape(theta',1,grid),reshape(theta',1,grid));
minD = abs(angdiff(Y,X));% .* deltatheta;
figure; imagesc(minD); axis square; 
colormap(copper); colorbar
title 'Dist_{min} between each E/I pairs ({\Delta}PD)'
xlabel 'x_{0^o,1^o,..,359^o}'
ylabel 'x_{0^o,1^o,..,359^o}'
set(gca,'YDir','Normal','YTick',[],'XTick',[])
%% Compute Weights
GsigmaDirE = kappaE * exp(-(minD.^2) / (2 * sigmadir^2));
GsigmaDirI = kappaI * exp(-(minD.^2) / (2 * sigmadir^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wee = rand(grid,grid);
Wee(Wee < GsigmaDirE) = 1;
Wee(Wee ~= 1) = 0;
Wee = Wee .* normrnd(Jee,0.25*Jee,[grid grid]);
Wee = max(Wee,0);
cWee = Jee * sum(GsigmaDirE,2);
aWee = sum(Wee,2);
sWee = cWee ./ aWee;
sWee(isinf(sWee)==1) = 0;
Wee = Wee .* repmat(sWee,1,grid);
% Wee = Jee .* GsigmaDir; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Wie = rand(grid,grid);
Wie(Wie < GsigmaDirE) = 1;
Wie(Wie ~= 1) = 0;
Wie = Wie .* normrnd(Jie,0.25*Jie,[grid grid]);
Wie = max(Wie,0);
cWie = Jie * sum(GsigmaDirE,2);
aWie = sum(Wie,2);
sWie = cWie ./ aWie;
sWie(isinf(sWie)==1) = 0;
Wie = Wie .* repmat(sWie,1,grid);
% Wie = Jie .* GsigmaDir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wii = rand(grid,grid);
Wii(Wii < GsigmaDirI) = 1;
Wii(Wii ~= 1) = 0;
Wii = Wii .* normrnd(Jii,0.25*Jii,[grid grid]);
Wii = max(Wii,0);
cWii = Jii * sum(GsigmaDirI,2);
aWii = sum(Wii,2);
sWii = cWii ./ aWii;
sWii(isinf(sWii)==1) = 0;
Wii = Wii .* repmat(sWii,1,grid);
% Wii = Jii .* GsigmaDir; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wei = rand(grid,grid);
Wei(Wei < GsigmaDirI) = 1;
Wei(Wei ~= 1) = 0;
Wei = Wei .* normrnd(Jei,0.25*Jei,[grid grid]);
Wei = max(Wei,0);
cWei = Jei * sum(GsigmaDirI,2);
aWei = sum(Wei,2);
sWei = cWei ./ aWei;
sWei(isinf(sWei)==1) = 0;
Wei = Wei .* repmat(sWei,1,grid);
% Wei = Jei .* GsigmaDir; 
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
xlabel 'x_{0^o,1^o,..,359^o}'
ylabel 'x_{0^o,1^o,..,359^o}'
suptitle 'Connection Weights'
%% h(theta)
hphi = ones(th_size,grid);
coh = 1; % stim coherence?
sigmaFF = deg2rad(60); % stimulus width
hphi = hphi .* (coh * exp(-(theta - deg2rad(180)).^2 / (2 * sigmaFF^2)));% + (1 - coh);
for i = 1:grid
    phi = theta(i); % stimulus direction phi
    hphi(i,:) = circshift(hphi(i,:),i-180);
end
figure;
imshow(hphi); axis square; colorbar
set(gca,'YDir','Normal')
title 'h(\theta)'
ylabel 'h(x_i)'
xlabel 'x_{0^o,1^o,..,359^o}'
%% Tau, k, & n
makeitRND = 'yes';
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
        uk = 0.04;
        k = normrnd(uk,0.05*uk, [2*grid 1]);
    case 'no'
        tau(mod(1:2*grid,2)==1) = tauE;
        tau(mod(1:2*grid,2)==0) = tauI;
        n = zeros(2*grid,1) + 2;
        k = zeros(2*grid,1) + 0.04;
end
return
%% Microstim
ms = 'on';
amp = -1;
switch ms
    case 'on'
        msExc = amp .* rand(grid,1);
%         hphi = meshgrid(msExc') + hphi; 
    case 'off'
        msExc = zeros(grid,1);
end
% Solve ODE 
c = 3;
stimpool = 90;
h = [hphi(stimpool,:)';hphi(stimpool,:)'];
msExc = [msExc;msExc];
% Solve ODE to obtain R
[t,r] = ode45(@(t,r) isn4(t,r,sparse(W),h,tau,k,c,n,msExc),TSpan,r0);
r = real(r);
Rei.fr = r;
Rei.t = t;
rf = r(end,:);

meanR = 'end';
switch meanR
    case 'end'
        FR = Rei(1).fr(end,:); % Response at the end of Evolution
    case 'all'
        FR = mean(Rei(1).fr);
end
rE = FR(1:360);%FR(mod(1:2*grid,2)==1);
rI = FR(361:end);%FR(mod(1:2*grid,2)==0);

%%
r = []; 
r(1,:,:) = rE;
r(2,:,:) = rI;
r = squeeze(mean(r));
ctl = 0:1:180;
R = cell(1,size(ctl,2));
hx_dir = ones(1,size(theta,2)) * theta(stimpool);
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
RmeanE = cellfun(@mean,RE);
RmeanI = cellfun(@mean,RI);

Rmean1m = Rmean;
%%
figure
plot(ctl,Rmean1m,'b','LineWidth',2); hold on
plot(ctl,Rmean0,'k','LineWidth',2)
plot(ctl,Rmean1,'r','LineWidth',2)
legend('-FB','Boring','+FB')
title 'mean E/I Response'
xlim([-5 185])
% ylim([0 1])
xlabel 'Deviation from MT PD^o'
ylabel 'Firing Rate'
%%
% Plot
figure
subplot(2,2,1)
plot(hphi(stimpool,:),'k','LineWidth',2)
legend 'h(\phi)'
xlabel 'E/I pair PD (deg)'
title(sprintf('Input Stim:  %d^o',stimpool))

pair = stimpool;
subplot(2,2,2)
plot(Rei.t,Rei.fr(:,pair-1),'b',Rei.t,Rei.fr(:,pair),'r','linewidth',2); hold on
legend('rE','rI')
box off
xlabel 'Time (ms)'
ylabel 'Firing rate'
title(sprintf('E/I pair@ %d^o',stimpool))

r = [];
% rE = FR(mod(1:2*grid,2)==1);
% rI = FR(mod(1:2*grid,2)==0);
r(1,:,:) = rE;
r(2,:,:) = rI;
r = squeeze(mean(r));
ctl = 0:1:180;
R = cell(1,size(ctl,2));
hx_dir = ones(1,size(theta,2)) * theta(stimpool);
EI_dev = floor(rad2deg(abs(angdiff(hx_dir,theta))));
for j = 1:size(ctl,2)
    R{j} = [R{j} r(EI_dev == ctl(j))];
end
Rmean = cellfun(@nanmean,R);

subplot(2,2,3)
plot(rE,'b','LineWidth',2); hold on
plot(rI,'r','LineWidth',2)
legend('rE','rI')
ylim([0 1])
xlabel 'E/I pair PD (deg)'
ylabel 'Firing rate'

% rE = FR(mod(1:2*grid,2)==1);
% rI = FR(mod(1:2*grid,2)==0);
ctl = 0:1:180;
RE = cell(1,size(ctl,2));
RI = RE;
hx_dir = ones(1,size(theta,2)) * theta(stimpool);
EI_dev = floor(rad2deg(abs(angdiff(hx_dir,theta))));
for j = 1:size(ctl,2)
    RE{j} = [RE{j} rE(EI_dev == ctl(j))];
    RI{j} = [RI{j} rI(EI_dev == ctl(j))];
end
RmeanE = cellfun(@mean,RE);
RmeanI = cellfun(@mean,RI);

subplot(2,2,4)
plot(ctl,Rmean,'k','LineWidth',2); hold on
plot(ctl,RmeanE,'b','LineWidth',1)
plot(ctl,RmeanI,'r','LineWidth',1)
legend('E/I','rE','rI')
% legend('rE','rI')
xlim([-5 185])
ylim([0 1])
xlabel 'Deviation from MT PD^o'
ylabel 'Firing Rate'
