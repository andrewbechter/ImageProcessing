%Process Image Data
% see readme for example code, and all the parameters/methods you have access to

%========================================%
% fitParameters format
% fitParameters(:,1) = amplitude
% fitParameters(:,2) = x centroid
% fitParameters(:,3) = x sigma
% fitParameters(:,4) = y centroid
% fitParameters(:,5) = y sigma
% fitParameters(:,6) = theta 
% fitParameters(:,7) = offset (i.e. background) 

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\ANDORV2\Aus4';
% load(filename)

filename ='P:\iLocater\NIC\ProcessedData\2017_07_07\SXFast2';
load(filename)
x = SXFast2.fitParameters(:,2) - SXFast2.Mean(1,2);
y = SXFast2.fitParameters(:,4) - SXFast2.Mean(1,4);
x = x*6.5*0.604;%(microns/pix)*(mas/micron)
y = y*6.5*0.604;
SXFast2 = calcPSD(SXFast2,y);
SXFast2.PSDPlot
% theta = atan(x./y);
% r = sqrt(x.^2+y.^2);
SXFast2 = calcPSD(SXFast2,y);
f = SXFast2.PSD(:,1);
bin = mean(diff(f));
SXFast2.PSD(:,2) = sqrt(SXFast2.PSD(:,2));
SXFast2.PSDPlot

figure
plot(SXFast2.time(:,1)*1e-6,x)
hold on
plot(SXFast2.time(:,1)*1e-6,x,'.')


grid on
xlabel('time (s)')
ylabel('Displacement (\mum)')
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 16;
hline = refline(0,std(x));
hline.Color = 'k';
hline.LineWidth = 2;
hline.LineStyle = '-.';

hline2 = refline(0,-std(x));
hline2.Color = 'k';
hline2.LineWidth = 2;
hline2.LineStyle = '-.';
