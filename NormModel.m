% clear all; clc

ctlmt = 0:45:315;
ctlmst = 0:45:315;
theta_i = deg2rad(ctlmst(randi(size(ctlmst,2))));
ncell = 20;

figure
for i = 1:ncell
    theta_j(i,1) = deg2rad(ctlmt(randi(size(ctlmt,2))));
    theta = abs(angdiff(theta_i,theta_j(i,1))) + deg2rad(90);
    
    params.bw = 2.5;
    params.pd = theta;
    params.fr = rand;
    params.spnt = -1;
    dir = deg2rad(0:45:315);
    

    D(i,:) = vonMises(params,dir,0);
    
    bw = [0.5,1,2];
    x0 = 20 * rand;
    y0 = 40 * rand - 20;
    
    if x0 <= 7 && (y0 >= -7 && y0 <= 7)
        r = bw(1);
    elseif (x0 >= 14 && x0 <= 20) && ((y0 >= -20 && y0 <= -14) || (y0 >= 14 && y0 <= 20))
        r = bw(3);
    else
        r = bw(2);
    end
    
    x = -20:0.5:20;
    y = x;
    [X,Y] = meshgrid(x,y);
    G(i,:,:) = Gauss2D(X,Y,x0,y0,r,1);
    hold on
end
line([-20 20],[0 0],'Color','k','LineWidth',1); hold on
line([0 0],[-20 20],'Color','k','LineWidth',1)
title 'Visual Field'
xlabel 'x-axis'
ylabel 'y-axis'

figure
e = 0.5;
for i = 1:ncell
    subplot(5,4,i)
    plot(rad2deg(dir),D(i,:),'k','LineWidth',2)
    xlim([min(rad2deg(dir))-e max(rad2deg(dir))+e])

end
xlabel 'Direction^{\circ}'
ylabel 'Firing Rate'

g = sum(sum(G,3),2);
d = mean(D,2);
w = sum(g .* d);


D = D ./ max(D,[],2);
D(D<0) = 0;

figure
for i = 1:20
    ch = conv2(D(i,:),squeeze(G(i,:,:)));
    subplot(5,4,i)
%     plot([0:45:315],D(i,:))
    mesh(ch)
%     pcolor(ch)
    w(i,1) = max(ch(:));
end
%%
R = [];
for i = 1:20
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
%%
a = [1;1;1;1;1];
w = [0.5;0.5;0.5;0.5;0.5];

wI = a .* w;
sigma = sum(wI);
R = a ./ sigma;
figure
plot(a,'k','LineWidth',2); hold on
plot(R,'r','LineWidth',2)








