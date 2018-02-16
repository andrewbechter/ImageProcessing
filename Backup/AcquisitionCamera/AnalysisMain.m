% Frquency Stuff

r = sqrt((I.Centroid(:,2)-I.Mean(1,2)+ I.Centroid(:,1)-I.Mean(1,2)).^2);

time_units = 1e-6; %seconds
dT = diff(I.Time)./diff(I.FrameNumbers).*time_units;
I.Frequency = 1./mean(dT);
[I]=Fourier_Transform (I,I.Centroid(:,2)-I.Mean(1,2));% Fourier transform of filtered or unfiltered data
[psor,lsor] = findpeaks(I.FFT(:,2),I.FFT(:,1),'SortStr','descend');
if isempty(psor)==1
    I.PeakFrq = 'None';
else
I.PeakFrq = lsor(1:5);
end
% Plotting
% I.PeakPlot(1)
% I.HistPeakPlot(1)
I.ScatterPlot
I.FFTPlot
I.PeakFrq