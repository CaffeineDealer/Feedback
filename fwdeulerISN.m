function [rout] = fwdeulerISN(Wee,Wei,node,h) %(tst,tend,Wee,Wei,Wie,Wii,node)

tst = 0;
tend = 400;
tsteps = 1; % 1ms fwd euler steps
N = (tend - tst) / tsteps; % # steps


re = zeros(size(Wee,1),1)+0.0001;
ri = re;


switch node
    case 'E'
        sigmaWR = (Wee .* re) - (Wei .* ri);
        sigmaWR = sum(sigmaWR(:));
        r = re;
        tau = 20;
        a = 1;
        hp = h;
        k = 1;
        n = 1;
        c = 1;
    case 'I'
        sigmaWR = sum((Wie * re) - (Wii .* ri));
        sigmaWR = sum(sigmaWR(:));
        a = 0;
        r = ri;
        tau = 10;
end

%%
% F = @(t,r) (-r + k.*(sigmaWR + c.*h + a.*hp).^n) / tau;

r(1) = 0.001;
rout = r;
tout = 0;
for t = tst:tsteps:tend
    sigmaWR = sigmaWR * .1;
    F = @(t,r) (-r + k.*(sigmaWR + c.*h + a.*hp).^n) / tau;
    r = r + h .* F(t,r);
    rout(:,t+1) = r;
end

figure
subplot(1,3,1)
plot(rout(:,end))
xlabel 'Network Position(degree)'
ylabel 'Firing Rate'
subplot(1,3,2)
plot(rout)
subplot(1,3,3)
plot(rout')
