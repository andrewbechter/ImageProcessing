% Image Main
tic

%% Fits or .fit files (Simulated files or SBIG)
% file = '10mm';
% ext = '.fit';
% path = 'P:\iLocater\GSFC_Fibers\NA2\1050nm\9_1\';
% filename = [path,file,ext];

%% Sifx(Andor)
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_15\BackIllumination\FlexureTest_Up\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_14\Pollux\Pollux_StackTest\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_13\Del02\Del02_Coupled\Spooled files.sifx';

%KxVir
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\KXVIR\KXVIR_sub\Spooled files.sifx';

%NuCom
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_sub_400\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NuCom\NuCom_full\Spooled files.sifx';

%Australis
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_153_sub\Spooled files.sifx';

%NIC
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_23_2\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_27_final\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\usb_floated_stability_06_28_3khz\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\2017_07_09\SX_center_1\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\2017_07_08\DX_Center\Spooled files.sifx';

%% Basler (.tiff or other standard image files. only specify the directory for Basler data)
% filename= 'P:\iLocater\NIC\Stability_NotFloated\Set1\';
% filename = 'P:\iLocater\NIC\Stability_Floated\Set2\'; 
% filename = 'P:\iLocater\QuadCell\Set1';
% filename = 'P:\Hardware_Calibration\WhiteLase\BeamProfile\1064_10x';
% filename = 'P:\iLocater\F41_Simulator\Image_data\Focus\Set10';
% filename = 'P:\iLocater\WhiteLase Coupler\Set1';
% filename = 'P:\iLocater\NA Spot Size test\Focus\Pair\Andrew\24.5mm';
% filename = 'P:\iLocater\NA Spot Size test\Focus\Pair\Ryan';
% filename = 'P:\iLocater\ImagingBoard\Set1';
% filename = 'P:\iLocater\QuadCell\Vis\Set23';
% filename = 'P:\iLocater\AcquisitionCamera\COTS\FiberChannel\Set16';
% filename = 'P:\iLocater\AcquisitionCamera\COTS\Andor\Set1';
filename = 'P:\iLocater\QuadCell\IR\Set109';

%% Simulated files
% filename = 'S:\SimulatorV9\TestScripts\AiryTest.mat';

%% Telescope sample files
% filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
% filename = 'Z:\2017_07_07\SX_vibration\Spooled files.sifx';

%% Create Image Object
I= Image;

%User Inputs
FitForOrientation=0; % set to 1 for rotated gaussin, 0 for non-rotated.
lowMem = 0; %Lowmem 0 to store frames, 1 to process without storage 
time_units = 1e-6; %seconds

%File reads
if length(strfind(filename,'.sifx')) == 1 %Andor sifx file
    [I.NoFrames] = Image.calc_num_frames(filename);
    mem_step = 1000;
    step = 1;
    I.step = mem_step;
    range = [1,10000];%I.NoFrames]; %range must start with 1 right now
%     range = [1,4500];
    frames = (range(1,1):step:range(1,2)); % frames to be processed using 1 index(ANDOR is 0 index)
    [frame_time,image_data,cuts,stored_frames] = Image.img_process_V2(filename,lowMem,FitForOrientation,frames,mem_step);
    I.cuts = cuts;
    I.Frame = double(image_data); I.Time = frame_time(:,1); I.stored_frames = stored_frames;%Assign properties to obj
    [I]= Process_to_object(I,frame_time); %Assign properties to obj using low mem format
    Test_frame =1;
    clear image_data frame_time
elseif length(strfind(filename,'.fits')) == 1 || length(strfind(filename,'.fit')) == 1 %Single fits file
    [I,frame] = SBIG_Frames(I,filename,FitForOrientation); %Basler directory
    [I]= Process_to_object(I,frame); %Assign properties to obj
    Test_frame =1;
elseif length(strfind(filename,'.mat')) == 1
    [I,frame] = Sim_Frames(I,filename,FitForOrientation);
    [I] = Process_to_object(I,frame);
    Test_frame = 1;
else
    [I,frame] = Basler_Frames(I,filename,FitForOrientation); %Basler directory
    time_units = 1; %seconds
    I.step =1;
    T = 25e-03; %data set time stamps (can't find any automated way to do this)
    L = size(frame,1); % Length of signal
    frame(:,1) = (0:L-1)*T; % overwrite time into the frame parameter variable 
    [I]= Process_to_object(I,frame); %Assign properties to obj
    I.NoFrames = length(I.FrameNumbers);
    Test_frame =1;
    clear L frame T 
end

% Plotting
I.PeakPlot(Test_frame)
toc

% BadRow removal
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

clear ii

if length (I.FrameNumbers)<50  %Not enough data for FFT and Histrograms
    return
end

% % Frquency Stuff is only needed with time series of frames
% dT = diff(I.Time)./diff(I.FrameNumbers).*time_units;
% % dT = 0.025;
% I.Frequency = 1./mean(dT);
% [I]=Fourier_Transform (I,I.Centroid(:,1));% Fourier transform of filtered or unfiltered data
% [psor,lsor] = findpeaks(I.FFT(:,2),I.FFT(:,1),'SortStr','descend');
% if isempty(psor)==1
%     I.PeakFrq = 'None';
% else
%     I.PeakFrq = lsor(1:5);
% end
I.HistPlot
% I.FFTPlot
I.Mean(1,[1,2])


