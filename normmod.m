clear all; clc; 
% close all
load('E:\MT_MST\Microstim\Model\mt.mat','TCNt');
Nt1 = TCNt;
clear TCNt;
load('E:\MT_MST\Microstim\Model\mst.mat','TCNt');
Nt2 = TCNt;
clear TCNt;

mt1 = [];
mt2 = [];
r1 = 1; r2 = 0;
for i = 1:size(Nt1,2)
    a = size(Nt2(i).MST,2);
    csize(i,1) = a;
    csize(i,2) = size(Nt1(i).MST.bMT,2);
    for j = 1:a
        mt1 = [mt1 Nt1(i).MST.bMT];
        r2 = r2 + csize(i,2);
        info(r1:r2,1) = i; % session #
        info(r1:r2,2) = j; % MST in session i
        info(r1:r2,3) = 1:csize(i,2); % MT in session i        
        r1 = r2 + 1;
    end
    mt2 = [mt2 Nt2(i).MST.bMT]; 
end
%% Arrange Firing Rates of MT
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
X1 = X1';
X2 = X2';
%%
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
%% cmass data import
load(['E:\MT_MST\Microstim\Model\' 'cmassMST.mat'])
cmMST(:,1) = Vx(:,2);
cmMST(:,2) = Vy(:,2);
Vx = []; Vy = [];
load(['E:\MT_MST\Microstim\Model\' 'cmassMT.mat'])
cmMT(:,1) = Vx(:,2);
cmMT(:,2) = Vy(:,2);
Vx = []; Vy = [];
% cmass in radians
cmMST = deg2rad(cmMST);
cmMT = deg2rad(cmMT);
%% MST & MT cmass
% MST cmass
for i = 1:22
    info(cnt(i,1):cnt(i,2),4) = cmMST(i,1); % x cmass MST
    info(cnt(i,1):cnt(i,2),5) = cmMST(i,2); % y cmass MST    
end

r1 = 1; r2 = 0;
for i = 1:6
    r2 = r2 + csize(i,2);
    for j = 1:csize(i,1)
        info(info(:,1)==i & info(:,2)==j,6) = cmMT(r1:r2,1); % x cmass MT
        info(info(:,1)==i & info(:,2)==j,7) = cmMT(r1:r2,2); % y cmass MT
    end
    r1 = r2 + 1;
end
%% MST & MT PD
load(['E:\MT_MST\Microstim\Model\' 'pd.mat'])
info(:,8) = pd(:,1); % PD MST
info(:,9) = pd(:,2); % PD MT
info(:,8) = deg2rad(indx2dir(info(:,8)));
info(:,9) = deg2rad(indx2dir(info(:,9)));  
%% Compute weights based on FC policy
load(['E:\MT_MST\Microstim\Model\' 'pdDiff.mat'])
dtheta = pdDiff(:,2) + deg2rad(90);
% dtheta = abs(angdiff(info(:,9),info(:,8))) + deg2rad(90);
dcmass = info(:,4:7);

params.bw = 2.5;
params.pd = dtheta*-1;
params.fr = X1(:,1);
params.spnt = 0;
dir = 0;
D = vonMises(params,dir,0);
bw = deg2rad(10);
G = Gauss2D(dcmass(:,1),dcmass(:,2),dcmass(:,3),dcmass(:,4),bw,0);
w = D .* G;

eps = 0.1;
figure
subplot(3,1,1)
plot(D,'.','MarkerSize',10,'Color','b')
title(sprintf('D_{VonMises}| bw = %0.1f, Mineault et al.,2012',params.bw))
ylim([min(D)-eps max(D)+eps])
subplot(3,1,2)
plot(G,'.','MarkerSize',10,'Color','r')
title(sprintf('G_{2DGaussian}| bw = %.1f, Reynolds & Heeger 2009',bw))
ylim([min(G)-eps max(G)+eps])
subplot(3,1,3)
plot(w,'.','MarkerSize',10,'Color','k')
title 'w = D.G| Minoldger et al., 2010.5'
ylim([min(w)-eps max(w)+eps])
%% Sigma w computation & assignment
for i = 1:size(csize,1)
    a = reshape(w(info(:,1)==i,1),csize(i,2),csize(i,1));
    a = sum(a,2);
    w(info(:,1)==i,2) = repmat(a,csize(i,1),1);    
    info(info(:,1)==i,10) = repmat(a,csize(i,1),1);
    a = [];
end
%% Normalization Model
% figure
pm.M = 1;
pm.sigma = 0.5;
for i = 1:size(csize,1)
    for j = 1:csize(i,1)
        q = find(info(:,1) == i & info(:,2) == j);
        pm.w = w(q,2);
        res(j,:,:) = fbnorm(X1(q,:),X2(q,:),pm);
    end
    X(i,:) = mean(X1(q,:),1);
    R(i,:) = mean(squeeze(mean(res)),1);
    res = [];
%     subplot(3,2,i)
%     plot(mean(X1(q,:),1),'k','LineWidth',2)
%     hold on
%     plot(R(i,:),'r','LineWidth',2)
end
%% Just some fucking graphs
mX = mean(X,1);
mR = mean(R,1);

errX = (std(X,[],1)) / (sqrt(size(X,1)));
errR = (std(R,[],1)) / (sqrt(size(R,1)));

figure
col(1) = 'k';
col(2) = 'r';
dir = 0:45:180;
plot(dir,mean(X,1),col(1),'LineWidth',2)
hold on
plot(dir,mean(R,1),col(2),'LineWidth',2)
xlabel 'Dev. from MT PD'
ylabel 'Norm. Firing Rate'
legend({'Observed','Modeled'})


shade(dir,mX+errX,sprintf('--%s',col(1)),dir,mX-errX,sprintf('--%s',col(1)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(1))
shade(dir,mR+errR,sprintf(':%s',col(2)),dir,mR-errR,sprintf(':%s',col(2)),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col(2))
%%
load(['E:\MT_MST\Microstim\Model\mt_data.mat'])

mN = mean(N,2);
mS = mean(S,2);
obs = mN - mS;
obt = mX - mR;

ob = N - S;
mo = X - R;

errob = (std(ob,[],2)) / (sqrt(size(ob,2)));
errmo = (std(mo,[],1)) / (sqrt(size(mo,1)));

ymin = min([min(obs) min(obt)]);

col1 = 'b';
col2 = 	'r';
figure
plot(dir,obs,'LineWidth',2,'Color',col1); hold on
plot(dir,obt,'LineWidth',2,'Color',col2)
ylim([ymin-0.1 1])
xlabel 'Dev. from MT PD'
ylabel 'Norm. diff Firing Rate'
legend({'Observed','Modeled'})

shade(dir,obs+errob,sprintf('--%s',col1),dir,obs-errob,sprintf('--%s',col1),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col1)
shade(dir,obt+errmo,sprintf(':%s',col2),dir,obt-errmo,sprintf(':%s',col2),'FillType',[1 2;1 2],'FillAlpha',.1,'FillColor',col2)

