%MainProcessingScript
clear
tic
%function descriptions:
%name: Andor
%purpose: grabs a single frame and time stamp
%inputs: ## filaname: filename including directory and .sifx
%        ## memoryStep: interval to store frames (storing frames uses up memory so be careful. The total frames stored should be ~100 max)
%        ## frameEnd: frame number to end on (if specified to be [] or more than the number frames
%           in spool, then all frames are processed).  
%        ## type: 'fast' or 'full' processing. fast uses approximate parameters, full uses 2D Gaussian fit
%        
%notes: only the filename needs to be specified. defaults are : no frames
%stored, all frames processed using 'full'
%outputs: Andor object with fit parameters, time stamps and stored frames
%if specified. 

%name: analyzeAndorData
%purpose: calculates useful things using time stamps and fit parameters
%inputs: ## Andor object
%notes: this method modifies exisiting Andor objects by populating or
%overwriting useful things for analysis 
%outputs: ## Andor object
tic
%%------------------
% Processing Optics
%%------------------

%NIC plate scale  = 0.604 ''/mm

memStep = 1000; % save interval

startFile = 1; % start frame

endFile = 0; % end frame

opts.backgroundfile = false; %backgroundfile

opts.trim = 300; %intial trim

opts.fineGrid = 256; %fine trim (min 256 for demonstrator b/s strehl computation)

opts.fitGrid = 25; %fine trim

opts.threshold = 1; % min counts threshold

%%------------------
% Optical parameters
%%------------------

opts.lambda = 1064e-9; % wavelegnth in meters

opts.focal_length = 250e-3; % effective focal length

opts.d = 6.068e-3; % beam diameter

opts.pix = 6.5e-6; %pixel size

opts.npup = 256*10; %number of samples across the pupil

opts.alpha = 0; %blocking fraction

opts.instSR = true; % compute instantaneous SR

opts.srgrid = 64; %computation grid set for SR copmutation

opts.fitPSF = true; % compute fitted PSF

q = opts.focal_length*opts.lambda/(opts.d*opts.pix);

tic

%=============================================%
%SX
% filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
% SX = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% SX = analyzeAndorData(SX,opts); %calculates standard analysis
% save Results/New/SX SX
% clear SX filename

filename = 'P:\iLocater\NIC\2017_07_08\DX_fast2\Spooled files.sifx';
DX2 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
DX2 = analyzeAndorData(DX2,opts); %calculates standard analysis
save Results/New/DX2 DX2
clear DX2 filename

toc