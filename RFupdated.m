function [Vx Vy diffV] = RFupdated(S,N)





posInf = struct;
for i = 1:size(S,2)
    posInf.fmaxposS(i,:) = max(S(i).firing); % Max firing rate between 8 directions @ each 9 position
    posInf.fmaxposN(i,:) = max(N(i).firing);
end
fmaxS = vertcat(posInf.fmaxposS);
fmaxN = vertcat(posInf.fmaxposN);
frmax(:,1) = max(fmaxS');
frmax(:,2) = max(fmaxN');
%%
%%%% Center of Mass %%%%
sepration = 17;
diam = 20;
x = [-17;0;17;-17;0;17;-17;0;17];
y = [17;17;17;0;0;0;-17;-17;-17];
for i = 1:size(fmaxS,1)
    [Vx(i,1) Vy(i,1)] = cmass(fmaxS(i,:)',x,y);
    [Vx(i,2) Vy(i,2)] = cmass(fmaxN(i,:)',x,y);
end
diffV(:,1) = Vx(:,1) - Vx(:,2); % x-axis
diffV(:,2) = Vy(:,1) - Vy(:,2); % y-axis
%%%% RF center of mass %%%%
% figure
% scatter(Vx(:,1),Vy(:,1),'FaceColor','b'); hold on
% scatter(Vx(:,2),Vy(:,2),'FaceColor','r')
% line([0 0],[min(y) max(y)],'LineWidth',0.5,'Color','k')
% line([min(x) max(Vx(:))],[0 0],'LineWidth',0.5,'Color','k')
% xlim([min(x) max(Vx(:))])
% line([min(x) max(x)],[0 0],'LineWidth',0.5,'Color','k')
% xlim([min(x) max(x)])
% ylim([min(y) max(y)])
% xlabel 'Position^o'
% ylabel 'Position^o'
