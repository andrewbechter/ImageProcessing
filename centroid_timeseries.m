%centroid time series figure
clear all

load('S:\ImageProcessing\Results\New\Aus5.mat')
x =  (Aus5.fitParameters(:,2)-Aus5.Mean(2))*2.33;
y =  Aus5.fitParameters(:,4);
t =  Aus5.time(:,1)*1e-6; % convert time into seconds
[px,Sx] = polyfit(t,x,1);
[py,Sy] = polyfit(t,y,1);
fit = polyval(px,t);

figure
plot(t,x)
hold on
plot(t,fit,'linewidth',1.5,'linestyle','-.')
xlim([0,10])
xlabel('time (s)')
ylabel('Centroid position (mas)')
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1;