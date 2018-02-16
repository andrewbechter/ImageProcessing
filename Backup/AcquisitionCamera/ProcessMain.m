% Image Main
clear
tic

%Read file (this should be setup such that the filenames for ANDOR and
%Power meter are paired together, (maybe vibrations and temperature too)

% Fits or .fit files (Simulated files or SBIG)
% filename = 'Y:\NIC_Photon.fits';

% Sifx (Andor)
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_15\BackIllumination\FlexureTest_Up\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_sub_400\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_full\Spooled files.sifx';

% filename = 'P:\iLocater\NIC\usb_floated_stability_06_23_2\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_27_final\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_28_3khz\Spooled files.sifx';

% Basler (.tiff or other standard image files. only specify the directory for Basler data)
% filename= 'P:\iLocater\NIC\Stability_NotFloated\Set1\';
% filename = 'P:\iLocater\NIC\Stability_Floated\Set2\'; 

% filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
% filename = 'Z:\2017_07_07\SX_vibration\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';

%Create Image Object
I= Image;

%User Inputs
FitForOrientation=1; % set to 1 for rotated gaussin, 0 for non-rotated.
lowMem = 1; %Lowmem 0 to store frames, 1 to process without storage 
time_units = 1e-6; %seconds

%File reads
if length(strfind(filename,'.sifx')) == 1 %Andor sifx file
    [I.NoFrames] = Image.calc_num_frames(filename);
    step = 1;
    range = [1,I.NoFrames]; %range must start with 1 right now
    range = [1,30000];
    frames = (range(1,1):step:range(1,2)); % frames to be processed using 1 index(ANDOR is 0 index)
    [frame_time,image_data] = Image.img_process(filename,lowMem,FitForOrientation,frames);
    I.Frame = double(image_data); I.Time = frame_time(:,1); %Assign properties to obj
    [I]= Process_to_object(I,frame_time); %Assign properties to obj using low mem format
    Test_frame =1;
    clear image_data frame_time
elseif length(strfind(filename,'.fits')) == 1 || length(strfind(filename,'.fit')) == 1 %Single fits file
    [I.Frame]=fitsread(filename);
    
    
else
    [I] = Basler_Frames(I,filename); %Basler directory
    time_units = 1; %seconds
    T = 0.001; %data set time stamps (can't find any automated way to do this)
    L = length(I.Frame); % Length of signal
    I.Time(:,1) = (0:L-1)*T;
end

%BadRow removal
[badInd] = find(I.Flag == 1);
I.Centroid(badInd,:)=[];
I.ApproxCentroid(badInd,:)=[];
I.Time(badInd,:) = [];
I.Flag(badInd,:) = [];
I.Frame(:,:,badInd)=[];

% Statistics
for ii = 1:size(I.Centroid,2)
I.Mean(1,ii) = mean(I.Centroid(:,ii)); % Mean value of Centroid parameters
I.RMS(1,ii) = std(I.Centroid(:,ii)); % RMS value of Centroid parameters
end

for ii = 1:2
I.Range(:,ii) = max(I.Centroid(:,ii))-min(I.Centroid(:,ii)); % Range
I.Residuals(:,ii) = abs(I.Centroid(:,ii))- min(I.Centroid(:,ii)); % Residuals
end

% Frquency Stuff
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
I.PeakPlot(Test_frame)
% I.HistPeakPlot(Test_frame)
I.ScatterPlot
I.FFTPlot
toc