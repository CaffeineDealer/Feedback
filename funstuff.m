clear all; clc
x = -2:0.02:2;
y = x;
f = 2.*x.*(cos(x.^2+y.^2));
g = 2.*y.*(cos(x.^2+y.^2));
[X,Y] = meshgrid(f,g);
Z = sin(X.^2 + Y.^2);
h = figure;
imagesc(Z)
colormap('copper')
% axis square
set(gca,'YDir','Normal')
h.Color = 'w';
box off
axis off
%%
clear all; clc
% webcamlist
cam = webcam;
preview(cam)
tic
c = 1;
while toc < 5
    a(:,:,:,c) = snapshot(cam);
    c = c + 1;
end
closePreview(cam)

figure
for i = 1:100
    subplot(10,10,i)
    colormap('gray')
    imagesc(a(:,:,1,i))
end
close all