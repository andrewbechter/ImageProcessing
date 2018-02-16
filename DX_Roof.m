load('S:\ImageProcessing\Reduced\DX_08_07_17\DX_iLocater.mat')
iloc_DX = I.Frame;
load('S:\ImageProcessing\Reduced\DX_08_07_17\DX_LMIRCam.mat')
LMIR_DX = I.Frame;

total = LMIR_DX + iloc_DX;

pscale = 1; % arcseconds/mm
pix = 6.5e-3; %mm
X = (1:size(LMIR_DX,2))*pix*pscale;
Y = (1:size(LMIR_DX,1))*pix*pscale;
X = X - mean(X);
Y = Y-mean(Y);

figure
imagesc(X,Y,total)
axis equal
axis image
title('DX Roof Mirror')
set(gca,'YDir','normal')
xlabel('\Delta mm')
ylabel('\Delta mm')