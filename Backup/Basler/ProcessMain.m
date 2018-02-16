% Image Main
clear all
tic
%Read file (this should be setup such that the filenames for ANDOR and
%Power meter are paired together, (maybe vibrations and temperature too)

%% Fits or .fit files (Simulated files or SBIG)
% filename = 'Y:\NIC_Photon.fits';
% filename = 'P:\iLocater\GSFC_Fibers\NA\1064_mid.fit';

% file = '10mm';
% ext = '.fit';
% path = 'P:\iLocater\GSFC_Fibers\NA\HighRes3\';
% filename = [path,file,ext];

% Sifx(Andor)
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_15\BackIllumination\FlexureTest_Up\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_sub_400\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_full\Spooled files.sifx';

% filename = 'P:\iLocater\NIC\usb_floated_stability_06_23_2\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_27_final\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_28_3khz\Spooled files.sifx';

%% Basler (.tiff or other standard image files. only specify the directory for Basler data)
% filename= 'P:\iLocater\NIC\Stability_NotFloated\Set1\';
% filename = 'P:\iLocater\NIC\Stability_Floated\Set2\'; 
% filename = 'P:\iLocater\QuadCell\Set1';
% filename = 'P:\Hardware_Calibration\WhiteLase\BeamProfile\1064_10x';

%% Telescope sample files
% filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
% filename = 'Z:\2017_07_07\SX_vibration\Spooled files.sifx';

%% Create Image Object
I= Image;

%User Inputs
FitForOrientation=0; % set to 1 for rotated gaussin, 0 for non-rotated.
lowMem = 1; %Lowmem 0 to store frames, 1 to process without storage 
time_units = 1e-6; %seconds

%File reads
if length(strfind(filename,'.sifx')) == 1 %Andor sifx file
    [I.NoFrames] = Image.calc_num_frames(filename);
    mem_step = 100;
    step = 1;
    I.step = mem_step;
    range = [1,I.NoFrames]; %range must start with 1 right now
%     range = [1,1000];
    frames = (range(1,1):step:range(1,2)); % frames to be processed using 1 index(ANDOR is 0 index)
    [frame_time,image_data,cuts] = Image.img_process(filename,lowMem,FitForOrientation,frames,mem_step);
    I.cuts = cuts;
    I.Frame = double(image_data); I.Time = frame_time(:,1); %Assign properties to obj
    [I]= Process_to_object(I,frame_time); %Assign properties to obj using low mem format
    Test_frame =1;
    clear image_data frame_time
elseif length(strfind(filename,'.fits')) == 1 || length(strfind(filename,'.fit')) == 1 %Single fits file
    [I,frame] = SBIG_Frames(I,filename,FitForOrientation); %Basler directory
    [I]= Process_to_object(I,frame); %Assign properties to obj
    Test_frame =1;
else
    [I,frame] = Basler_Frames(I,filename,FitForOrientation); %Basler directory
    time_units = 1; %seconds
    T = 0.001; %data set time stamps (can't find any automated way to do this)
    L = size(frame,1); % Length of signal
    frame(:,1) = (0:L-1)*T; % overwrite time into the frame parameter variable 
    [I]= Process_to_object(I,frame); %Assign properties to obj
    I.NoFrames = length(I.FrameNumbers);
    Test_frame =1;
    clear L frame T 
end

%BadRow removal
% [badInd] = find(I.Flag == 1);
% [badInd] = find(I.FitResidual >0.1);
% I.Centroid(badInd,:)=[];
% I.ApproxCentroid(badInd,:)=[];
% I.Time(badInd,:) = [];
% I.Flag(badInd,:) = [];
% I.Frame(:,:,badInd)=[];

% Statistics
for ii = 1:size(I.Centroid,2)
I.Mean(1,ii) = mean(I.Centroid(:,ii)); % Mean value of Centroid parameters
I.RMS(1,ii) = std(I.Centroid(:,ii)); % RMS value of Centroid parameters
end

for ii = 1:2
I.Range(:,ii) = max(I.Centroid(:,ii))-min(I.Centroid(:,ii)); % Range
I.Residuals(:,ii) = abs(I.Centroid(:,ii))- min(I.Centroid(:,ii)); % Residuals
end
clear ii

% Plotting
I.PeakPlot(Test_frame)
toc

% Frquency Stuff is only needed with time series of frames
if length (I.FrameNumbers)<50   
    return
else
    dT = diff(I.Time)./diff(I.FrameNumbers).*time_units;
    I.Frequency = 1./mean(dT);
    [I]=Fourier_Transform (I,I.Centroid(:,2));% Fourier transform of filtered or unfiltered data
    [psor,lsor] = findpeaks(I.FFT(:,2),I.FFT(:,1),'SortStr','descend');
    if isempty(psor)==1
        I.PeakFrq = 'None';
    else
        I.PeakFrq = lsor(1:5);
    end
    I.HistPLot
    I.FFTPlot
end


