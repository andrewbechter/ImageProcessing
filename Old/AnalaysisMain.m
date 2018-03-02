% clear

% load('S:\ImageProcessing\Reduced\Australis\Australis6.mat')
% load('S:\ImageProcessing\Reduced\KxVir\KxVir.mat')
load('S:\ImageProcessing\Reduced\NuCom\NuCom.mat')


time_units = 1e-6;

% Frquency Stuff is only needed with time series of frames
dT = diff(I.Time)./diff(I.FrameNumbers).*time_units;
I.Frequency = 1./mean(dT);
I.FFT = [];
R = (sqrt((I.Centroid(:,1)-I.Mean(1,1))+(I.Centroid(:,2)-I.Mean(1,2))).^2);
[I]=Fourier_Transform (I,R*2.33);% Fourier transform of filtered or unfiltered data
[psor,lsor] = findpeaks(I.FFT(:,2),I.FFT(:,1),'SortStr','descend');

if isempty(psor)==1
    I.PeakFrq = 'None';
else
    I.PeakFrq = lsor(1:5);
end

% Analysis
Test_frame = 1;

% Plotting
I.PeakPlot(Test_frame)

% BadRow removal
% [badInd] = find(I.FitResidual >0.2);
[badInd1] = find(I.Centroid(:,1) > I.Mean(1)+4*I.RMS(1));
[badInd2] = find(I.Centroid(:,1) < I.Mean(1)-4*I.RMS(1));
[badInd3] = find(I.Centroid(:,2) > I.Mean(2)+4*I.RMS(2));
[badInd4] = find(I.Centroid(:,2) < I.Mean(2)-4*I.RMS(2));

badInd = [badInd1;badInd2;badInd3;badInd4];
% assert(length(badInd)<0.1*length(I.Centroid(:,1)),'too may indicies removed');

I.Centroid(badInd,:)=[];
I.ApproxCentroid(badInd,:)=[];
I.Time(badInd,:) = [];
I.FrameNumbers(badInd,:) = [];
I.Flag(badInd,:) = [];

% Statistics
for ii = 1:size(I.Centroid,2)
I.Mean(1,ii) = mean(I.Centroid(:,ii)); % Mean value of Centroid parameters
I.RMS(1,ii) = std(I.Centroid(:,ii)); % RMS value of Centroid parameters
end

I.Residuals = [];
I.Range = [];
for ii = 1:2
I.Range(:,ii) = max(I.Centroid(:,ii))-min(I.Centroid(:,ii)); % Range
I.Residuals(:,ii) = abs(I.Centroid(:,ii))- min(I.Centroid(:,ii)); % Residuals
end

I.AngHistPlot
I.FFTPlot




