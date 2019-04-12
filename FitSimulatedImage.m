%Fit simulated image
load('/Volumes/Software/PsfAberrations/PSF2.mat')

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

gridx = mean(abs(diff(FPgridx(1,:))));
gridy = mean(abs(diff(FPgridy(:,1))));
[x] = Image.subPixelPeakV3(frame,x0,FitForOrientation,xdata);

x(3) = x(3)*2.355;
x(5) = x(5)*2.355;
fprintf('Amp: %.5f\n x0: %.2f\n FWHMx: %.2f\n y0: %.2f\n FWHMy: %.2f\n theta: %.2f\n offset: %.2f\n',x) 