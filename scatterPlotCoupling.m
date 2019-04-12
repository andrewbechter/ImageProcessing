%Scatter plot with coupling efficiency
clear 
close 

load('/Volumes/Software/ImageProcessing/Results/New/NuComSub400.mat')

load('/Volumes/Software/Simulator/TiptiltStudy/tt_v_fc.mat')

f = fit(data(:,1),data(:,2),'gauss1');

x = (min(data(:,1)):max(data(:,1)));
[X,Y] = meshgrid(x,x');

exponent = (X.^2 + Y.^2)./(f.c1.^2);
Z = (exp(-exponent));

%surface(X,Y,Z)
%NuComSub400.masPlot

v = [0.9,0.5];

NuComSub400.masPlot
pause
grid on
[C,h] = contour(X,Y,Z,v,'color','k','linestyle','--','linewidth',2);
clabel(C,h,'FontSize',15)



% 
% r1 = 5 ;% mas radius equivilant to 0.9 coupling 
% [x1,y1] = circle(0,0,r1);
% 
% r2 = 10;% mas radius equivilant to 0.6 coupling 
% [x2,y2] = circle(0,0,r2);
% 
% r3 = 20;% mas radius equivilant to 0.3 coupling
% [x3,y3] = circle(0,0,r3);
% 
% hold on
% h = plot(x1, y1, 'linewidth', 1.5, 'linestyle','--','color','k' );
% 
% function [xunit,yunit] = circle(x,y,r)
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% end
