%AnalyzeAndor
clear all
tic
%function descriptions:
%name: AnalyzeAndor
%purpose: load Andor files and re-run analysis 

%-----------
%Load Files 
%-----------

load('S:\ImageProcessing\Results\New\NuComSub300.mat')
load('S:\ImageProcessing\Results\New\NuComSub400.mat')
load('S:\ImageProcessing\Results\New\Aus1.mat')
load('S:\ImageProcessing\Results\New\Aus2.mat')
load('S:\ImageProcessing\Results\New\Aus3.mat')
load('S:\ImageProcessing\Results\New\Aus4.mat')
load('S:\ImageProcessing\Results\New\Aus5.mat')
load('S:\ImageProcessing\Results\New\Aus6.mat')
load('S:\ImageProcessing\Results\New\Aus153.mat')
load('S:\ImageProcessing\Results\New\KXVIRSub.mat')
load('S:\ImageProcessing\Results\New\NSV.mat')
load('S:\ImageProcessing\Results\New\DX2.mat')
load('S:\ImageProcessing\Results\New\SX2.mat')


%----------------
%Analysis Options
%----------------

opts.backgroundfile = false; %backgroundfile

opts.trim = 300; %intial trim

opts.fineGrid = 256; %fine trim (min 256 for demonstrator b/s strehl computation)

opts.fitGrid = 25; %fine trim

opts.threshold = 1; % min counts threshold

opts.lambda = 800e-9; % wavelegnth in meters

opts.focal_length = 750e-3; % effective focal length

opts.d = 13.33e-3; % beam diameter

opts.pix = 6.5e-6; %pixel size

opts.npup = 256*10; %number of samples across the pupil

opts.alpha = 0.11; %blocking fraction

opts.instSR = true; % compute instantaneous SR

opts.srgrid = 256; %computation grid set for SR copmutation

opts.fitPSF = true; % compute fitted PSF

%----------------
%Analyze On-Sky
%----------------

Aus1 = analyzeAndorData(Aus1,opts); %calculates standard analysis
d(1,1) = Aus1.drift.*2.33.*60;
d(1,2) = Aus1.drift.*Aus1.time(end,1)*1e-6;
d(1,3) = Aus1.time(end,1)*1e-6;
close all

Aus2 = analyzeAndorData(Aus2,opts); %calculates standard analysis

d(2,1) = Aus2.drift.*2.33.*60;
d(2,2) = Aus2.drift.*Aus2.time(end,1)*1e-6;
d(2,3) = Aus2.time(end,1)*1e-6;

close all

Aus3 = analyzeAndorData(Aus3,opts); 
d(3,1) = Aus3.drift.*2.33.*60;
d(3,2) = Aus3.drift.*Aus3.time(end,1)*1e-6;
d(3,3) = Aus3.time(end,1)*1e-6;
close all

Aus4 = analyzeAndorData(Aus4,opts); 
d(4,1) = Aus4.drift.*2.33.*60;
d(4,2) = Aus4.drift.*Aus4.time(end,1)*1e-6;
d(4,3) = Aus4.time(end,1)*1e-6;
close all

Aus5 = analyzeAndorData(Aus5,opts);
d(5,1) = Aus5.drift.*2.33.*60;
d(5,2) = Aus5.drift.*Aus5.time(end,1)*1e-6;
d(5,3) = Aus5.time(end,1)*1e-6;
close all

Aus6 = analyzeAndorData(Aus6,opts); 
d(6,1) = Aus6.drift.*2.33.*60;
d(6,2) = Aus6.drift.*Aus6.time(end,1)*1e-6;
d(6,3) = Aus6.time(end,1)*1e-6;
close all

Aus153 = analyzeAndorData(Aus153,opts); 
d(7,1) = Aus153.drift.*2.33.*60;
d(7,2) = Aus153.drift.*Aus153.time(end,1)*1e-6;
d(7,3) = Aus153.time(end,1)*1e-6;
close all

NuComSub300 = analyzeAndorData(NuComSub300,opts); 
d(8,1) = NuComSub300.drift.*2.33.*60;
d(8,2) = NuComSub300.drift.*NuComSub300.time(end,1)*1e-6;
d(8,3) = NuComSub300.time(end,1)*1e-6;
close all

NuComSub400 = analyzeAndorData(NuComSub400,opts);
d(9,1) = NuComSub400.drift.*2.33.*60;
d(9,2) = NuComSub400.drift.*NuComSub400.time(end,1)*1e-6;
d(9,3) = NuComSub400.time(end,1)*1e-6;
close all

% NSV = analyzeAndorData(NSV,opts); 
% NSV.drift.*2.33.*60;

KXVIRSub = analyzeAndorData(KXVIRSub,opts); 
d(10,1) = KXVIRSub.drift.*2.33.*60;
d(10,2) = KXVIRSub.drift.*KXVIRSub.time(end,1)*1e-6;
d(10,3) = KXVIRSub.time(end,1)*1e-6;
% d(10) = 1.53;
close all

%----------------
%Closed Dome
%----------------

DX2 = analyzeAndorData(DX2,opts); 
d(11,1) = DX2.drift.*0.604.*6.5*60; %(mas/micron) *1 pix in microns
d(11,2) = DX2.drift.*DX2.time(end,1)*1e-6;
d(11,3) = DX2.time(end,1)*1e-6;
close all

SX2 = analyzeAndorData(SX2,opts); 
d(12,1)= SX2.drift.*0.604.*6.5*60;
d(12,2) = SX2.drift.*SX2.time(end,1)*1e-6;
d(12,3) = SX2.time(end,1)*1e-6;
close all

toc