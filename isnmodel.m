clear all; clc
%% Input
grid = 2; % grid size
%% MST & MT PD
load(['D:\MT_MST\Microstim\Model\' 'pdDiff.mat'])
mt = pdDiff(1:8,2);
%% The map of PD of the E/I pairs
% maybe the PD difference between each pair of MT?
% grid = size(mt,1); % E/I pairs
theta = zeros(grid,grid);
mt = reshape(mt,grid,grid,2);
theta = abs(angdiff(mt(:,:,1),mt(:,:,2)));

figure;
imagesc(theta); axis square; colorbar
%% Parameters
deltax = min(theta(theta>0)); % angular distance sepration of neurons to be shortest distance between deltaPD = min(theta)
% deltax = representation of visual space / map size;

%parameters
kappaE = 0.1;
kappaI = 0.5;
Jee = 0.044;
Jie = 0.042;
Jei = 0.023;
Jii = 0.018;

sigmaEE = 12 * deltax;
sigmaIE = sigmaEE;
sigmaEI = 4 * deltax;
sigmaII = sigmaEI;
sigmadir = deg2rad(64); % What is this?
sigmaRF = deltax;
%% Weights
% Let Wab(x,x') be the synaptic weight form the cell of type b, position
% x', preferred orientation theta', to the cell of type a, position x,
% preferred orientation theta. Nonzero connections are sparse and chosen
% randomly, with probability

[x,y] = meshgrid(1:grid,1:grid);
X = [reshape(x,1,grid^2) ; reshape(y,1,grid^2)]; 
xyD = distmat(X');
xyD = xyD * deltax;

[X,Y] = meshgrid(reshape(theta',1,grid^2),reshape(theta',1,grid^2));
minD = abs(X-Y);
% minD(minD>180) = 360 - minD(minD>180);

pWE = kappaE * exp(-xyD.^2/(2*sigmaEE^2)).*exp(-(minD).^2/(2*sigmadir^2));

Wee = rand(grid^2,grid^2);
Wee(Wee<pWE) = 1;
Wee(Wee~=1) = 0;
Wee = Wee.*normrnd(Jee,0.25*Jee,[grid^2 grid^2]);
Wee = max(Wee,0);
figure; imagesc(Wee); axis square; colorbar

cWee = Jee*sum(pWE,2);
aWee = sum(Wee,2);
sWee = cWee./aWee;
sWee(isinf(sWee)==1)=0;
Wee = Wee.*repmat(sWee,1,grid^2);

Wie = rand(grid^2,grid^2);
Wie(Wie<pWE) = 1;
Wie(Wie~=1) = 0;
Wie = Wie.*normrnd(Jie,0.25*Jie,[grid^2 grid^2]);
Wie = max(Wie,0);
cWie = Jie*sum(pWE,2);
aWie = sum(Wie,2);
sWie = cWie./aWie;
sWie(isinf(sWie)==1) = 0;
Wie = Wie.*repmat(sWie,1,grid^2);

pWI = kappaI*exp(-xyD.^2/(2*sigmaEI^2)).*exp(-(minD).^2/(2*sigmadir^2));

Wii = rand(grid^2,grid^2);
Wii(Wii<pWI) = 1;
Wii(Wii~=1) = 0;
Wii = Wii.*normrnd(Jii,0.25*Jii,[grid^2 grid^2]);
Wii = max(Wii,0);
figure; imagesc(Wii); axis square; colorbar

cWii = Jii*sum(pWI,2);
aWii = sum(Wii,2);
sWii = cWii./aWii;
sWii(isinf(sWii)==1) = 0;
Wii = Wii.*repmat(sWii,1,grid^2);

Wei = rand(grid^2,grid^2);
Wei(Wei<pWI) = 1;
Wei(Wei~=1) = 0;
Wei = Wei.*normrnd(Jei,0.25*Jei,[grid^2 grid^2]);
Wei = max(Wei,0);
cWei = Jei*sum(pWI,2);
aWei = sum(Wei,2);
sWei = cWei./aWei;
sWei(isinf(sWei)==1) = 0;
Wei = Wei.*repmat(sWei,1,grid^2);

% Where a nonzero connection exists, Wab(x,x') is chosen randomly from a
% Gaussian distirbution with mean Jab and standard deviation 0.25Jab;
% weights below zero are set to zero.

% For each cell, the set of recurrent synaptic weights of type b (E
% or I it receives are then scaled so that all cells of a given
% type a (E or I) receive the same total inhibitory and the same
% total excitatory synaptic weight form the newtwork, equal to Jab
% times the mean number of connections received under
% p(Wab(x,x'~=0).

W=[Wee -Wei; Wie -Wii];

%%
%input h(x)=sl(x) (this describes a slice through a diameter of the 2D circularly symmetric stimulus)

h = zeros(grid,grid);
cen = round(grid/2);
l = 20;
l = l*deltax;
[X,Y] = meshgrid(1:1:grid, 1:1:grid);

coh = 1;
sigmaFF = 2*17.623;
phi = theta(cen,cen);
itheta = reshape(theta',1,grid^2);
g = coh*exp(-(itheta-phi).^2/(2*sigmaFF^2))+(1-coh);

h = (1./(1+exp(-(sqrt((X-cen).^2+(Y-cen).^2)*deltax+l/2)/sigmaRF))).*(1-1./(1+exp(-(sqrt((X-cen).^2+(Y-cen).^2)*deltax-l/2)/sigmaRF)));
h = reshape(h,1,grid^2).*g;

figure
imagesc(reshape(h,grid,grid)');axis square; colorbar

h = [h';h'];

%contrast
c = 1;
%%

TSpan = [0 400];
r0 = zeros(2*grid^2,1);
% When I is set to 0
%r0=[0;0];
% In all cases, all drawn from Gaussian distributions, with standard
% deviation 0.05 times the mean
utauE = 20;
utauI = 10;
tau = [normrnd(utauE,0.05*utauE, [grid^2 1]);normrnd(utauI,0.05*utauI, [grid^2 1])];
% Nonlinearity, imaginery parts ignored?
unE = 2.0;
unI = 2.2;
n=[normrnd(unE,0.05*unE, [grid^2 1]);normrnd(unI,0.05*unI, [grid^2 1])];
uk=0.012;
k=normrnd(uk,0.05*uk, [2*grid^2 1]);

%%
[t,r]=ode45(@(t,r) isn4(t,r,sparse(W),h,tau,k,c,n),TSpan,r0);
%%

figure
handle=plot(t,r(:,grid^2-1),'b',t,r(:,grid^2),'r','linewidth',2);
legend('rE','rI')
box off
xlabel('Time (ms)','fontsize',12,'fontweight','b')
ylabel('Firing rate','fontsize',12,'fontweight','b')

%excitatory and inhibitory neuron
FR=r(end,:);
%FR=FR/Max(FR);
rE=FR(mod(1:2*grid^2,2)==1);
rI=FR(mod(1:2*grid^2,2)==0);

figure
rE=reshape(rE,grid,grid)';
imagesc(real(rE));axis square; colorbar

rI=reshape(rI,grid,grid)';
figure
imagesc(real(rI));axis square; colorbar
