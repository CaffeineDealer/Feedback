clear all; clc; %close all
%% This is just to make a grid and then superpose an orientation map on the grid 
grid = 360;    % grid size
th_size = 360;
theta = zeros(1,th_size);
theta = deg2rad(linspace(0,359,th_size));
%% Compute Direction preference for each pair
deltatheta = deg2rad(1);
[X,Y] = meshgrid(reshape(theta',1,grid),reshape(theta',1,grid));
minD = abs(angdiff(Y,X));% .* deltatheta;
% figure; imagesc(minD); axis square; 
% colormap(copper); colorbar
% title 'Dist_{min} between each E/I pairs ({\Delta}PD)'
% xlabel 'x_{0^o,1^o,..,359^o}'
% ylabel 'x_{0^o,1^o,..,359^o}'
% set(gca,'YDir','Normal','YTick',[],'XTick',[])
%% parameters
Jee = 0.044; % unstable 20.044
Jie = 0.042; 
Jei = 0.023; 
Jii = 0.018; % unstable 20.018
sigmadir = deg2rad(64);
kappaE = 0.1;
kappaI = 0.2;
if Jei*Jie > Jee*Jii
    disp('Stable')
else
    disp('Unstable')
end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
phi = 180; % stim dir
hphi = ones(1,grid);
coh = .5; % stim coherence? it was 1
sigmaFF = deg2rad(60); % stimulus width - it was 30
hphi = hphi .* (coh * exp(-(theta - deg2rad(180)).^2 / (2 * sigmaFF^2)));
hphi = circshift(hphi,phi-180);
% FB
FBe = ones(grid,1);
FBi = FBe;
% Feedback Parameters
sigmaFBe = deg2rad(90); % width it was 60
sigmaFBi = deg2rad(40); % width it was 60
muE = .3; % amplitude it was .3
muI = .3; % amplitude it was .3
% Excitatory Feedback
FBe = FBe' .* (muE * exp(-(abs(angdiff(theta,ones(1,360)*deg2rad(phi)))).^2 / (2 * sigmaFBe^2)));
% Inhibitory Feedback
FBi1 = FBi' .* (muI * exp(-(abs(angdiff(theta,ones(1,360)*(deg2rad(phi)+(pi/2))))).^2 / (2 * sigmaFBi^2)));
FBi2 = FBi' .* (muI * exp(-(abs(angdiff(theta,ones(1,360)*(deg2rad(phi)-(pi/2))))).^2 / (2 * sigmaFBi^2)));
FBi = FBi1 + FBi2;
FB = [FBe';FBi'];
% Plot section
figure;
plot(hphi,'k','LineWidth',2); hold on
plot(FBe,'r','LineWidth',2)
plot(FBi,'b','LineWidth',2)
legend({'h(x)','E','I'})
xlim([0 360])
ylim([0 1.1])
xlabel 'x_{0^o,1^o,..,359^o}'
st = trapz(hphi); 
fbe = trapz(FBe); 
fbi = trapz(FBi); 
disp(sprintf('stim = %.2f',st))
disp(sprintf('fbe = %.2f',fbe))
disp(sprintf('fbi = %.2f',fbi))
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
        n = ones(2*grid,1)*2; % n can range between 1.6<n<2.4
        k = ones(2*grid,1)*.04; % .04
end
%% Solve ODE
FB = [-FBe,-FBi];
a = [.5;1];
b = [{'FB_E + FB_I'},{'0_E + FB_I'},{'FB_E + 0_I'}];
figure
plot(FBe,'r','LineWidth',2); hold on
plot(FBi,'b','LineWidth',2)
legend({'E','I'},'box','off')
xlim([0 360])
xlabel('E/I Pair Pd^{o}','FontSize',15,'FontName','High Tower Text')
ylabel('Amplitude_{p.u.}','FontSize',15,'FontName','High Tower Text')
set(gca,'Box','off')
set(gca,'XTick',[0 45 90 135 180 225 270 315 360],'XTickLabel',[0 45 90 135 180 225 270 315 359])
%%
figure
for i = 1:2
    [t,r] = ode45(@(t,r) isn2E(t,r,sparse(W),[hphi';hphi'],tau,1,k,n,a(i),FB'),TSpan,r0);
    Fr(i,:) = real(r(end,:));
%     Fr(i,:) = mean(real(r),1);
    [Rmean(i,:) ctl E(i,:) I(i,:)] = fr2dev(Fr(i,:),'end',theta,phi,grid);
end
Rmean = Rmean ./ max(Rmean(1,:));
plot(ctl,Rmean(1,:),'k','LineWidth',2); hold on
plot(ctl,Rmean(2,:),'Color',[0 .5 0],'LineWidth',2)
xlim([0 180])
set(gca,'XTick',[0 45 90 135 180],'XTickLabel',[0 45 90 135 180])
set(gca,'Box','off')
legend({'Control','MS'},'box','off')
xlabel('\DeltaPD_{MT}^{o}','FontSize',15,'FontName','High Tower Text')
ylabel('Firing Rate_{p.u.}^{E/I Pairs}','FontSize',15,'FontName','High Tower Text')
set(gcf,'color','w')

% E = E ./ max(E,[],2);
% figure
% plot(ctl,E(1,:),'k','LineWidth',2); hold on
% plot(ctl,E(2,:),'Color',[0 .5 0],'LineWidth',2)
return
%%
clear all; clc
n = 1.6;
Jee = 0.044;
Jii = 0.018;
for i = 1
    params.efb = 1;
    params.d = 0.1;
    params.n = n;
    params.Jee = Jee;
    params.Jii = Jii;
%     [e(i,:,:),~,~] = ssn_run(params);
%     E(i,1) = e(i,1,1);
    [~,~,e(i,:,:)] = ssn_run(params);
    n = n - 0;
    Jee = Jee + .3;
    Jii = Jii + .3;
end
%%
figure
for i = 1
%     subplot(1,5,i)
    R = squeeze(e(i,:,:));
    R = R ./ max(R(1,:));
    
    plot(R(1,:),'k','LineWidth',1.5); hold on
    hold on
    plot(R(2,:),'Color',[0 .5 0],'LineWidth',1.5)
    xlim([0 180])
    set(gca,'XTick',[0 45 90 135 180],'XTickLabel',[0 45 90 135 180])
    set(gca,'Box','off')
    if i == 1
        legend({'Control','MS'},'box','off')
        xlabel('\DeltaPD_{MT}^{o}','FontSize',15,'FontName','High Tower Text')
        ylabel('Firing Rate_{p.u.}^{E/I Pairs}','FontSize',15,'FontName','High Tower Text')
        set(gcf,'color','w')
    end
end
return
%%
figure;
semilogy(E,'k','LineWidth',2) 
hold on; box off
semilogy(e(:,2,1),'g','LineWidth',2)
set(gca,'Box','off')
legend({'Control','MS'},'box','off')
xlabel('Jee.Jii: Stable \rightarrow Unstable ','FontSize',15,'FontName','High Tower Text')
ylabel('Excitatory Cell Firing Rate_{\DeltaPD = 0^{o}}','FontSize',15,'FontName','High Tower Text')
set(gcf,'color','w')


