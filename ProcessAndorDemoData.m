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

memStep = 1000; % save interval

startFile = 1; % start frame

endFile = 1; % end frame

% opts.delta = 9; % sets the frame size for fitting (2*delta)

opts.backgroundfile = false; %backgroundfile

opts.trim = 300; %intial trim

opts.fineGrid = 256; %fine trim (min 256 for demonstrator b/s strehl computation)

opts.fitGrid = 25; %fine trim

opts.threshold = 1; % min counts threshold

%%------------------
% Optical parameters
%%------------------

opts.lambda = 800e-9; % wavelegnth in meters

opts.focal_length = 750e-3; % effective focal length

opts.d = 13.33e-3; % beam diameter

opts.pix = 6.5e-6; %pixel size

opts.npup = 256*10; %number of samples across the pupil

opts.alpha = 0.11; %blocking fraction

opts.instSR = true; % compute instantaneous SR

opts.srgrid = 256; %computation grid set for SR copmutation

opts.fitPSF = true; % compute fitted PSF

%=============================================%
%AUSTRALIS

filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx';
Aus1 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus1 = analyzeAndorData(Aus1,opts); %calculates standard analysis
% save Results/New/Aus1 Aus1
% clear Aus1 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_2\Spooled files.sifx';
% Aus2 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus2 = analyzeAndorData(Aus2,opts); %calculates standard analysis
% save Aus2 Aus2
% clear Aus2 filename 
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_3\Spooled files.sifx';
% Aus3 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus3 = analyzeAndorData(Aus3,opts); %calculates standard analysis
% save Aus3 Aus3
% clear Aus3 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_4\Spooled files.sifx';
% Aus4 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus4 = analyzeAndorData(Aus4,opts); %calculates standard analysis
% save Aus4 Aus4
% clear Aus4 filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_5\Spooled files.sifx';
% Aus5 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus5 = analyzeAndorData(Aus5,opts); %calculates standard analysis
% save Aus5 Aus5
% clear Aus5 filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';
% Aus6 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus6 = analyzeAndorData(Aus6,opts); %calculates standard analysis
%save Results/New/Aus6 Aus6
% clear Aus6 filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_153_sub\Spooled files.sifx';
% Aus153 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus153 = analyzeAndorData(Aus153,opts); %calculates standard analysis
% save Results/New/Aus153 Aus153
% clear Aus153 filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_153_full\Spooled files.sifx';
% Aus_full = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus_full = analyzeAndorData(Aus_full,opts); %calculates standard analysis
% save Aus Aus
% clear Aus filename

%=============================================%
%72LEO
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\72Leo\72Leo_full\Spooled files.sifx';
% Leo72Full = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Leo72Full = analyzeAndorData(Leo72Full,opts); %calculates standard analysis
% save Leo72Full Leo72Full
% clear Leo72Full filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\72Leo\72Leo_sub\Spooled files.sifx';
% Leo72Sub = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Leo72Sub = analyzeAndorData(Leo72Sub,opts); %calculates standard analysis
% save Leo72Sub Leo72Sub
% clear Leo72Sub filename

%=============================================%
%KXVIR
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\KXVIR\KXVIR_full\Spooled files.sifx';
% test = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% test = analyzeAndorData(test,opts); %calculates standard analysis
% save KXVIRFull KXVIRFull
% clear KXVIRFull filename
% 

% endFile = 3278; % data bad after this frame
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\KXVIR\KXVIR_sub\Spooled files.sifx';
% KXVIRSub = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% KXVIRSub = analyzeAndorData(KXVIRSub,opts); %calculates standard analysis
% save Results/New/KXVIRSub KXVIRSub
% clear KXVIRSub filename

%=============================================%
%NuCom
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NUCOM\NuCom_Full\Spooled files.sifx';
% NuComFull = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% NuComFull = analyzeAndorData(NuComFull,opts); %calculates standard analysis
% save NuComFull NuComFull
% clear NuComFull filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NUCOM\nuCom_Sub_400\Spooled files.sifx';
% NuComSub400 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% NuComSub400 = analyzeAndorData(NuComSub400,opts); %calculates standard analysis
% save Results/New/NuComSub400 NuComSub400
% clear NuComSub400 filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\NUCOM\nuCom_Sub\Spooled files.sifx';
% NuComSub300 = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% NuComSub300 = analyzeAndorData(NuComSub300,opts); %calculates standard analysis
% save Results/New/NuComSub300 NuComSub300
% clear NuComSub filename

%=============================================%
%NSV
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_14\NSV19434\NSV19434\Spooled files.sifx';
% NSV = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% NSV = analyzeAndorData(NSV,opts); %calculates standard analysis
% save Results/New/NSV NSV
% clear NuComSub filename

toc