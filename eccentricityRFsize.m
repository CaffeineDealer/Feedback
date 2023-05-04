% v3 RF mapping
names={'slu023b','slu022e', 'slu023b','slu047d', 'slu050b', 'slu055c', 'slu062b', 'slu060c'};
path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
count=1; countv3=0; countmt=0;
count2=1;
probesize=10;
centroidV3=[]; xEccV3=[]; yEccV3=[];
sizeV3=[]; xDiaV3=[];yDiaV3=[];
centroidMT=[]; xEccMT=[]; yEccMT=[];
sizeMT=[]; xDiaMT=[]; yDiaMT=[]; dayV3=[]; dayMT=[];
errorV3=[]; errorMT=[]; ipsV3=[]; ipsMT=[]; ipsV32=[]; ipsMT2=[];
OLV3=[];OLMT=[];
pv3=[]; pmt=[];
ipsIdxRF=[1 2 7 8 13 14 19 20]
timeWin=(0.05*Fs:0.4*Fs);
for j=1:length(names)
    load(['C:\research\V3 things\V3 categorized\',names{1,j}(1:end-1),'_V3categ2.mat']);
    v3categ=sortrows(v3categ2);
    V3units=v3categ((v3categ(:,3)<=4),1:2);
    MTunits=v3categ(v3categ(:,3)==5,1:2);
    for ci=1:size(V3units,1)+size(MTunits,1)
        if ci<=size(V3units,1)
            ch=V3units(ci,1);
            u=V3units(ci,2);
            countv3=countv3+1;
        else
            ch=MTunits(ci-size(V3units,1),1);
            u=MTunits(ci-size(V3units,1),2);
            countmt=countmt+1;
        end
        spktrain=load([path,names{1,j},num2str(ch),num2str(u),'spktrain.mat']);
        spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
        baseline=squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
        allstimfir=squeeze(sum(spktrain.spktrain,1))*Fs/size(spktrain.spktrain,1);
        
        [h,p] = ttest(baseline(:),allstimfir(:));
        keepCriteria=(p<=0.05)&&max(allstimfir(:))>3;
        if keepCriteria
            firing=squeeze(mean(sum(spktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin),5));
            bestdir=max(firing,[],1);
            for jjj=1:4
                RF(jjj,:)=bestdir((jjj-1)*6+1:jjj*6);
            end
            x = -1*(6-1)*10/2:10:(6-1)*10/2;
            y = -1*[-1*(4-1)*10/2:10:(4-1)*10/2];
            %% thresh
            thresh=mean(baseline(:));
            RFth=(RF>thresh);
            [peak,I] = max(RF(:));
            RFth=(RF>0.6*peak);
            Size=sqrt(sum(RFth(:))*10*10);
            [Px,Py]=ind2sub([4,6],I);
            centroid=sqrt((x(Py))^2+(y(Px))^2);
            ips=x(Py);
            %% interpolation
            RFinterp= interp2(RF,5);
            [m,n]=size(RFinterp);
            [I,J]=ndgrid(-m/2:m/(m-1):m/2,-n/2:n/(n-1):n/2);
            correctionX=40/m;
            correctionY=60/n;
            [fitresult, zfit, fiterr, zerr, resnorm, rr] = fmgaussfit(I,J,RFinterp);
            centroid2=sqrt((fitresult(5)*correctionX)^2+(fitresult(6)*correctionY)^2);%*probesize;
            xEcc=fitresult(5)*correctionX;
            yEcc=fitresult(6)*correctionY;
            sigma2=(abs(fitresult(3)*correctionX/sqrt(2))+abs(fitresult(4)*correctionY/sqrt(2)))/2; %average diameter heuer 1999
            Size2=sigma2;%*probesize;
            xDia=abs(fitresult(3)*correctionX/sqrt(2));
            yDia=abs(fitresult(4)*correctionY/sqrt(2));
            ips2=(fitresult(5)*correctionX);
            overlap2=-(fitresult(5)*correctionX)-abs(fitresult(3)*correctionX/sqrt(2));
            %%
            firing1=sum(spktrain.spktrain(timeWin,:,:,ipsIdxRF,:,:,:),1)*Fs/length(timeWin);
            firingbl=sum(spktrainbl.spktrain_bl(timeWin,:,:,ipsIdxRF,:,:,:),1)*Fs/length(timeWin);
            for l=1:length(ipsIdxRF)
                x2=firing1(:,:,:,l,:,:,:);
                y2=firingbl(:,:,:,l,:,:,:);
                [h1,p2(l),c1,stats]=ttest(x2(:),y2(:))
            end
            p3=min(p2);
            significant=p3<(0.05/8)
            
            %%
            if ci<=size(V3units,1)
                centroidV3=[centroidV3,centroid2];
                xEccV3=[xEccV3, xEcc];
                yEccV3=[yEccV3,yEcc];
                xDiaV3=[xDiaV3, xDia];
                yDiaV3=[yDiaV3,yDia];
                %sizeV3=[sizeV3,Size-14.14]; %2*HALF OF THE SQUARE DIAGONAL * 10 DEG
                sizeV3=[sizeV3,Size2];
                errorV3=[errorV3, rr];
                ipsV3=[ipsV3, ips]
                ipsV32=[ipsV32, ips2]
                OLV3=[OLV3, overlap2]
                pv3=[pv3 significant]
                dayV3=[string(dayV3) names{1,j}]
            else
                centroidMT=[centroidMT,centroid2];
                xEccMT=[xEccMT, xEcc];
                yEccMT=[yEccMT, yEcc];
                xDiaMT=[xDiaMT, xDia];
                yDiaMT=[yDiaMT,yDia];
                %sizeMT=[sizeMT,Size-14.14];
                sizeMT=[sizeMT,Size2];
                errorMT=[errorMT, rr];
                ipsMT=[ipsMT, ips]
                ipsMT2=[ipsMT2, ips2]
                OLMT=[OLMT, overlap2]
                pmt=[pmt significant]
                dayMT=[string(dayMT) names{1,j}]
            end
        end
    end
end
%%
ignore=(sizeV3>=60|sizeV3<5)|(errorV3>=2000)|(centroidV3<5);
sizeV3r=sizeV3(~ignore);
centroidV3r=centroidV3(~ignore);
xEccV3r=xEccV3(~ignore);
yEccV3r=yEccV3(~ignore);
xDiaV3r=xDiaV3(~ignore);
yDiaV3r=yDiaV3(~ignore);
errorV3r=errorV3(~ignore);
dayV3r=dayV3(~ignore)
pv3r=pv3(~ignore)
sum(pv3r)/length(pv3r)

ignore2=(sizeMT>=60|sizeMT<5)|(errorMT>=2000)|(centroidMT<5);
sizeMTr=sizeMT(~ignore2);
centroidMTr=centroidMT(~ignore2);
xEccMTr=xEccMT(~ignore2);
yEccMTr=yEccMT(~ignore2);
xDiaMTr=xDiaMT(~ignore2);
yDiaMTr=yDiaMT(~ignore2);
errorMTr=errorMT(~ignore2);
dayMTr=dayMT(~ignore2)
pmtr=pmt(~ignore2)
sum(pmtr)/length(pmtr)

bV3=[ones(length(centroidV3r),1) centroidV3r']\sizeV3r';
PMT=polyfit(centroidMTr,sizeMTr,1);
PV3=polyfit(centroidV3r,sizeV3r,1);
bMT=[ones(length(centroidMTr),1) centroidMTr']\sizeMTr';
%%
figure;
subplot (1,2,1)
scatter(centroidMTr,sizeMTr,'*')
hold on
scatter(centroidV3r,sizeV3r,'o','MarkerFaceColor',[0.87,0.49,0])
hold on
plot(centroidMTr,[ones(length(centroidMTr),1) centroidMTr']*bMT,'LineWidth',...
    1.5,'Color',[0,0.45,0.74])
hold on
plot(centroidV3r,[ones(length(centroidV3r),1) centroidV3r']*bV3,'LineWidth',...
    1.5,'Color',[0.87,0.49,0])
xlim([0 max([centroidV3r,centroidMTr])+2])
ylim([0 max([sizeV3r,sizeMTr,centroidV3r,centroidMTr])])
xticks([0, 10, 20, 30])
yticks([0, 10, 20, 30, 40])
legend('MT','V3A')%,'fit2 V3','fit2 MT')
xlabel('Eccentricity (deg.)')
ylabel('Size (deg.)')
subplot(1,2,2)
day='slu023b'%'slu023b'
idx1=find(dayV3r==day)
idx2=find(dayMTr==day)
for i=idx1%length(yEccV3r)
    DrawCircle(xEccV3r(i),yEccV3r(i),xDiaV3r(i)/2, yDiaV3r(i)/2,[0.87,0.49,0])
end
for k=idx2%length(yEccMTr)
    DrawCircle(xEccMTr(k),yEccMTr(k),xDiaMTr(k)/2,yDiaMTr(k)/2,[0,0.45,0.74])
end
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlim([-60 60])
ylim([-40 40])
xlabel('x position (deg.)')
ylabel('y position (deg.)')
%%
figure
scatter3(centroidV3r,sizeV3r,errorV3r)
%Isline()
eccent=[centroidV3r,centroidMTr];
SIZE=[sizeV3r,sizeMTr];
region=[zeros(length(centroidV3r),1);ones(length(centroidMTr),1)];
aoctool(eccent,SIZE,region)
[h,p] = ttest2(sizeV3r',sizeMTr')

OLThresh=[-20 -10 0 10 20 30 40 50 60];
for di=1:length(OLThresh)
    thresh=OLThresh(di);
    nmt(di)=length(OLMT(-OLMT>=thresh))/length(OLMT);
    nv3(di)=length(OLV3(-OLV3>=thresh))/length(OLV3);
end
[h,p]=ttest2(OLMT',OLV3')
% % Two-sample Kolmogorov-Smirnov test
[h,p] = kstest2(OLMT,OLV3)
figure
plot(OLThresh,nmt,'LineWidth',2)
hold on
plot(OLThresh,nv3,'LineWidth',2)
legend('MT','V3A')
xlabel('Ipsilateral Overlap (Deg)')
ylabel('Percentile')
ylim([0 1])


