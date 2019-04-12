% ProcessBasler
%clear 

%%------------------
% Processing Options
%%------------------
fs = 1; % sampling frequency

memStep = 1; % save interval

startFile = 1; % start frame

endFile = 100; % end frame

opts.trim = 1000; %initial trim

opts.srgrid = 64; %computation grid set for SR copmutation

opts.fineGrid = 64; %fine trim (min 256 for demonstrator b/c strehl computation)

opts.fitGrid = 5; %fit window (should be about the PSF size)

opts.threshold = 0; % min counts threshold

%%------------------
% Optical parameters
%%------------------

opts.lambda = 1064e-9; % wavelegnth in meters

opts.npup = 1000; % number of samples across the pupil

opts.alpha = 0; % blocking fraction

opts.instSR = true; % compute instantaneous SR

opts.fitPSF = true; % compute fitted PSF

%%------------------
% Fiber Channel
%%------------------

opts.focal_length = 22e-3; %22e-3 is nominal, %L2 = 20e-3; %L3 = 257e-3;% optimal = 25.7e-3;

opts.d = 6.068e-3; % beam diameter

opts.pix = 1.67e-6; % pixel size

q = opts.focal_length*opts.lambda/(opts.d*opts.pix);

%%------------------
% Simulator
%%------------------
% load('CalibrationFiles/14umBasler/Set2_sim_back')
% opts.focal_length = 228.6000e-3; %;
% opts.d = 5.55e-3; % beam diameter
% opts.pix = 1.67e-6; %4.8e-6;
 
opts.q = opts.focal_length*opts.lambda/(opts.d*opts.pix);

%%------------------
% Data frames
%%------------------
%filename = 'P:\iLocater\QuadCell\IR\Set225\';
%filename = 'P:\iLocater\WhiteLaseCoupler\Set13\';
%filename = '/Volumes/Projects/iLocater/AcquisitionCamera/COTS/FiberChannel/Set17/';
%filename = 'P:\iLocater\/AcquisitionCamera/COTS/FiberChannel/Set17/';
%filename = 'P:/iLocater/AcquisitionCamera/Final/Fiber/L1Set2/';
%filename = 'P:/iLocater/AcquisitionCamera/Final/Simulator/Set1/';

%%-----------------------------------
% Fiber channel 
% L1 + L2 (#2) reference with background
%%--------------------------------------
% filename = 'P:/iLocater/AcquisitionCamera/Final/Fiber/Set5/';
filename = 'P:/iLocater/AcquisitionCamera/Final/FiberChannel/Comb_Set6\';
load('CalibrationFiles/14umBasler/Comb_Set6_back')

%%----------------------------
% L1 reference with background
%%----------------------------
% filename = 'P:/iLocater/AcquisitionCamera/Final/Fiber/L1Set3/';

%%----------------------------
%L2 certification with background
%%----------------------------
% filename = 'P:/iLocater/AcquisitionCamera/Final/Fiber/FiberSet1_rings/';
% filename = 'P:/iLocater/AcquisitionCamera/Final/Fiber/L2Set8_1/';
% load('CalibrationFiles/14umBasler/L2Set8_back')

%%----------------------------
%Simulator reference 14um Basler with background
%%----------------------------
% filename = 'P:/iLocater/AcquisitionCamera/Final/Simulator/Set2/'; %

%%----------------------------
% ADC dispersion certification no background
%%----------------------------
% filename = 'P:\iLocater\ADC Performance\Single Prism Testing\prism2dispersion';


%%------------------
% Calibration Frame
%%------------------
% opts.backgroundfile = false; %backgroundfile

opts.backgroundfile = backgroundFrame;

tic
fiber = Basler(filename,opts,memStep,startFile,endFile,'full',fs); % reads in data and populates fundamental object properties
fiber = analyzeBaslerData(fiber,opts); %calculates standard analysis
%img.psfPlot
toc

