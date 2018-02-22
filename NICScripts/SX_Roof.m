load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_ilocater.mat')
iloc_SX = I.Frame;
load('S:\ImageProcessing\Reduced\SX_09_07_17\SX_LMIRCam.mat')
LMIR_SX = I.Frame;

total = LMIR_SX + iloc_SX;

pscale = 1; % arcseconds/mm
pix = 6.5e-3; %mm
X = (1:size(LMIR_SX,2))*pix*pscale;
Y = (1:size(LMIR_SX,1))*pix*pscale;
X = X - mean(X);
Y = Y-mean(Y);

figure
imagesc(X,Y,total)
axis equal
axis image
title('SX Roof Mirror')
set(gca,'YDir','normal')
xlabel('\Delta mm')
ylabel('\Delta mm')
