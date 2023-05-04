clear all; clc
%% This is just to make a grid and then superpose an orientation map on the grid 
grid=32;    % grid size

nm=30;  % Kaschube et al 2010 Supplementary equation 20, nm is wavevectors
kc=2*pi*2/32;
theta=zeros(grid,grid);

[X,Y]=meshgrid(1:1:grid, 1:1:grid);
for j=1:nm
    l=randi(2);
    if l==2
        l=-1;
    end
    phi=2*pi*rand(1);
    
    theta=theta+exp((l*(X*kc*cos((j-1)*pi/nm)+Y*kc*sin((j-1)*pi/nm))+phi).*1i);                     
end
theta=angle(theta)*180/pi+180;
figure; 
imagesc(theta); axis square; colorbar
title 'Map of Preferred Orientation of the E/I pair'
clear X Y
%% parameters
%The full map is taken to be 16/75*grid*16/75*grid deg
deltax = 16/75;

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
sigmadir = 64;
sigmaRF = deltax;
%% Compute the distance between grid intervals: regardless of E/I, the distance difference of each to grid interval is stored here
[x,y] = meshgrid(1:grid,1:grid);
X = [reshape(x,1,grid^2) ; reshape(y,1,grid^2)]; 
xyD = distmat(X');
xyD=xyD*deltax; %scale to deg
figure; imagesc(xyD); axis square
%% Compute Orientation preference for each pair
[X,Y] = meshgrid(reshape(theta',1,grid^2),reshape(theta',1,grid^2));
minD=abs(X-Y);
minD(minD>180)=360-minD(minD>180);
figure; imagesc(minD); axis square
%% Compute Weights
% Wee & Wie
pWE=kappaE*exp(-xyD.^2/(2*sigmaEE^2)).*exp(-(minD).^2/(2*sigmadir^2));

Wee=rand(grid^2,grid^2);
Wee(Wee<pWE)=1;
Wee(Wee~=1)=0;
Wee=Wee.*normrnd(Jee,0.25*Jee,[grid^2 grid^2]);
Wee=max(Wee,0);

cWee=Jee*sum(pWE,2);
aWee=sum(Wee,2);
sWee=cWee./aWee;
sWee(isinf(sWee)==1)=0;
Wee=Wee.*repmat(sWee,1,grid^2);
figure; imagesc(Wee); axis square

Wie=rand(grid^2,grid^2);
Wie(Wie<pWE)=1;
Wie(Wie~=1)=0;
Wie=Wie.*normrnd(Jie,0.25*Jie,[grid^2 grid^2]);
Wie=max(Wie,0);
cWie=Jie*sum(pWE,2);
aWie=sum(Wie,2);
sWie=cWie./aWie;
sWie(isinf(sWie)==1)=0;
Wie=Wie.*repmat(sWie,1,grid^2);
figure; imagesc(Wie); axis square

 % Wii & Wei
pWI=kappaI*exp(-xyD.^2/(2*sigmaEI^2)).*exp(-(minD).^2/(2*sigmadir^2));

Wii=rand(grid^2,grid^2);
Wii(Wii<pWI)=1;
Wii(Wii~=1)=0;
Wii=Wii.*normrnd(Jii,0.25*Jii,[grid^2 grid^2]);
Wii=max(Wii,0);

cWii=Jii*sum(pWI,2);
aWii=sum(Wii,2);
sWii=cWii./aWii;
sWii(isinf(sWii)==1)=0;
Wii=Wii.*repmat(sWii,1,grid^2);
figure; imagesc(Wii); axis square

Wei=rand(grid^2,grid^2);
Wei(Wei<pWI)=1;
Wei(Wei~=1)=0;
Wei=Wei.*normrnd(Jei,0.25*Jei,[grid^2 grid^2]);
Wei=max(Wei,0);
cWei=Jei*sum(pWI,2);
aWei=sum(Wei,2);
sWei=cWei./aWei;
sWei(isinf(sWei)==1)=0;
Wei=Wei.*repmat(sWei,1,grid^2);
figure; imagesc(Wei); axis square

W=[Wee -Wei; Wie -Wii];
figure; imagesc(W); axis square
%%
h=zeros(grid,grid);
cen=round(grid/2);
l=20;
l=l*deltax;
[X,Y]=meshgrid(1:1:grid, 1:1:grid);

coh=1;
sigmaFF=2*17.623;
phi=theta(cen,cen);
itheta=reshape(theta',1,grid^2);
g=coh*exp(-(itheta-phi).^2/(2*sigmaFF^2))+(1-coh);

h=(1./(1+exp(-(sqrt((X-cen).^2+(Y-cen).^2)*deltax+l/2)/sigmaRF))).*(1-1./(1+exp(-(sqrt((X-cen).^2+(Y-cen).^2)*deltax-l/2)/sigmaRF)));

h=reshape(h,1,grid^2).*g;

figure
imagesc(reshape(h,grid,grid)');axis square; colorbar
h=[h';h'];
%contrast
c=1;
%%
TSpan = [0 400];

% Set r(t=0)=I
r0 = zeros(2*grid^2,1);
% When I is set to 0
%r0=[0;0];

% In all cases, all drawn from Gaussian distributions, with standard
% deviation 0.05 times the mean
utauE = 20;
utauI = 10;
% tauE=normrnd(utauE,0.05*utauE, [grid^2 1]);
% tauI=normrnd(utauI,0.05*utauI, [grid^2 1]);
% tau=zeros(2*grid^2,1);
% tau(mod(1:2*grid^2,2)==1)=tauE;
% tau(mod(1:2*grid^2,2)==0)=tauI;
tau = [normrnd(utauE,0.05*utauE, [grid^2 1]);normrnd(utauI,0.05*utauI, [grid^2 1])];

% Nonlinearity, imaginery parts ignored?
unE = 2.0;
unI = 2.2;
% nE=normrnd(unE,0.05*unE, [grid^2 1]);
% nI=normrnd(unI,0.05*unI, [grid^2 1]);
% n=zeros(2*grid^2,1);
% n(mod(1:2*grid^2,2)==1)=nE;
% n(mod(1:2*grid^2,2)==0)=nI;
n = [normrnd(unE,0.05*unE, [grid^2 1]);normrnd(unI,0.05*unI, [grid^2 1])];

uk = 0.012;
k = normrnd(uk,0.05*uk, [2*grid^2 1]);

% n = zeros(2*grid^2,1);
% k = n;
% n = n + 2;
% k = k + 0.04;
%% Solve ODE to obtain R
[t,r]=ode45(@(t,r) isn4(t,r,sparse(W),h,tau,k,c,n),TSpan,r0);

figure
handle=plot(t,r(:,grid^2-1),'b',t,r(:,grid^2),'r','linewidth',2);
legend('rE','rI')
box off
xlabel('Time (ms)','fontsize',12,'fontweight','b')
ylabel('Firing rate','fontsize',12,'fontweight','b')
%% E & I neuron
FR=r(end,:);
rE=FR(mod(1:2*grid^2,2)==1);
rI=FR(mod(1:2*grid^2,2)==0);

figure
rE=reshape(rE,grid,grid)';
imagesc(real(rE));axis square; colorbar

rI=reshape(rI,grid,grid)';
figure
imagesc(real(rI));axis square; colorbar
