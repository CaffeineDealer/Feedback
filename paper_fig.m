
% sig = Xb{100}(35,:,4); %ytu316 329
sig = Xb{10}(10,:,2); %ytu329
t = linspace(-400,500,9000);
ulimit = max(sig);
llimit = min(sig);
h = figure;
plot(t,sig,'k'); hold on
line([0 0],[llimit ulimit],'Color',[0 .5 0],'LineWidth',2,'LineStyle','--')
ylim([llimit ulimit])
% box off
axis off
h.Color = 'w';
tstart = 5485;
tend = 7825;
% plot(t(tstart:tend),sig(tstart:tend),'Color',[0.6350, 0.0780, 0.1840],'LineWidth',1)
