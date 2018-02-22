clear 
close
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_LMIRCam.mat')
center = I.Frame;
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_corner1.mat')
corner1 = I.Frame;
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_corner2.mat')
corner2 = I.Frame;
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_corner3.mat')
corner3 = I.Frame;
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_corner4.mat')
corner4 = I.Frame;

total = center+corner1+corner2+corner3+corner4;
figure
pscale = 0.604; % arcseconds/mm
pix = 6.5e-3; %mm
X = (1:size(corner1,2))*pix*pscale;
Y = (1:size(corner1,1))*pix*pscale;
X = X - mean(X);
Y = Y-mean(Y);

imagesc(X,Y,total)
axis equal
axis image
title('SX off-axis FOV')
set(gca,'YDir','normal')
xlabel('\Delta arcseonds')
ylabel('\Delta arcseonds')
