%Main Processing Script for ANDOR Fine Guide Camera (FGC)

%function descriptions:

%%--------------------------
% name: Andor
%%--------------------------

% purpose: grabs a single frame and time stamp
% inputs: ## filaname: filename including directory and .sifx
%         ## opts: include processing options:
%         ## memoryStep: interval to store frames (storing frames uses up memory so be careful. The total frames stored should be ~100 max)
%         ## startFile: frame number to end on
%         ## endFile: frame number to end on (if specified to be [] or more than the number frames
%           in spool, then all frames are processed).  
%         ## type: 'fast' or 'full' processing. fast uses approximate parameters, full uses 2D Gaussian fit
%        
% notes: only the filename needs to be specified. defaults are : no frames
% stored, all frames processed using 'full'
% outputs: Andor object with fit parameters, time stamps and stored frames
% if specified. 

%%--------------------------
% name: analyzeAndorData
%%--------------------------

% purpose: calculates useful things using time stamps and fit parameters
% inputs: ## Andor object
% notes: this method modifies exisiting Andor objects by populating or
% overwriting useful things for analysis 
% outputs: ## Andor object

tic

%%------------------
% Processing Options
%%------------------

opts.memStep = 1; % save interval

opts. startFile = 1; % start frame

opts. endFile = 1; % end frame

opts.trim = 256; %intial trim

opts.fineGrid = 256; %fine trim (min 256 for demonstrator b/s strehl computation)

opts.fitGrid = 20; %fine trim

opts.threshold = 0; % min counts threshold

opts.type = 'full'; % processing type

%%------------------
% Optical parameters
%%------------------

opts.lambda = 1064e-9; % wavelegnth in meters

opts.focal_length = 250e-3; % effective focal length

opts.d = 6.068e-3; % beam diameter

opts.pix = 6.5e-6; %pixel size

opts.npup = 1000; %number of samples across the pupil

opts.alpha = 0; %blocking fraction

opts.instSR = true; % compute instantaneous SR

opts.srgrid = opts.fineGrid; %computation grid set for SR copmutation

opts.fitPSF = true; % compute fitted PSF

q = opts.focal_length*opts.lambda/(opts.d*opts.pix);

%%------------------
% Calibration Files
%%------------------

load('CalibrationFiles/Andor/190408_back')% load background file

opts.backgroundfile = backgroundFrame; %backgroundfile

%%------------------
% Data Files
%%------------------

% Final Acquisition camera
filename = 'P:\iLocater\AcquisitionCamera\Final\Andor\20190408\Strehl\Spooled files.sifx';

% PSD realtime
%filename = 'S:\ImageProcessing\MLA\PSD_realtime\Spooled files.sifx';

%%------------------
% Function Calls
%%------------------

FGC = Andor(filename,opts); % reads in data and populates fundamental object properties
FGC = analyzeAndorData(FGC,opts); %calculates standard analysis
% save Results/FGC

toc