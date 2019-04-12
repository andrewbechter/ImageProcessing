%Fit simulated image Andrew
% load('S:\ImageProcessing\SimulatedImages\test.mat')

APSF = abs(test4.PSF).^2;
dl = mean(abs(diff(test4.FPgridx(1,:))));

[X,Y] = meshgrid(1:size(APSF,1),1:size(APSF,2));
% Meshgrid steps over '0' which makes the frame 1 pixel larger than
% desired. Truncates to correct size. Offsets all frame values (i.e pixels)
% to be centered at location of PSF - represents the actual detector location
xdata(:,:,1) = X; %layer 1 is X search space
xdata(:,:,2) = Y; %layer 2 is Y search space

% xdata = cat(3,FPgridx,FPgridy) ;%meshgrid x and meshgrid y of grid data
FitForOrientation = 0;
x0 = [max(max(APSF)),size(APSF,2)/2,3/2.355,size(APSF,2)/2,3/2.355,0,0];%initial guess input format x0 = [Amp,xo,wx,yo,wy,theta,offset];
frame = APSF; %counts matrix


[x] = Image.subPixelPeakV3(frame,x0,FitForOrientation,xdata);


d = x;
d(2) = (x(2)-(size(APSF,1)/2+1))*dl*1e6;
d(3) = x(3)*dl*4*1e6;
d(4) = (x(4)-(size(APSF,2)/2+1))*dl*1e6;
d(5) = x(5)*dl*4*1e6;

fprintf('Amp: %.5f\n x0: %.2f\n FWHMx: %.2f\n y0: %.2f\n FWHMy: %.2f\n theta: %.2f\n offset: %.2f\n',d) 