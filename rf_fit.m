clear all; clc
% RF (fr_bestDir,pos)
RF = ones(3,3)*5;
RF(1,2) = 20;
% Interpolate RF
RFinterp = interp2(RF,5);
% RFinterp = RF;
figure
imagesc(RFinterp)
%%
[m,n] = size(RFinterp);
[I,J] = ndgrid(-m/2:m/(m-1):m/2,-n/2:n/(n-1):n/2);
[fitresult,zfit,fiter,zerr,resnorm,rr] = fmgaussfit(I,J,RFinterp);
correctionX = 40/m;
correctionY = 60/n;
centroid = sqrt((fitresult(5)*correctionX)^2+(fitresult(6)*correctionY)^2);%*probesize;
xEcc = fitresult(5)*correctionX;
yEcc = fitresult(6)*correctionY;
Sigma = (abs(fitresult(3)*correctionX/sqrt(2))+abs(fitresult(4)*correctionY/sqrt(2)))/2; %average diameter heuer 1999
xDia = abs(fitresult(3)*correctionX/sqrt(2));
yDia = abs(fitresult(4)*correctionY/sqrt(2));
ips2 = (fitresult(5)*correctionX);
overlap2 = -(fitresult(5)*correctionX)-abs(fitresult(3)*correctionX/sqrt(2));
%%
figure
scatter(centroid,Sigma,'o','MarkerFaceColor','k')
hold on
DrawCircle(xEcc,yEcc,xDia/2,yDia/2,'k')