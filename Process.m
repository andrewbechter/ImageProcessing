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

memStep = 10; % save interval

startFile = 1; % start frame

endFile = 0; % end frame

opts.backgroundfile = false; %backgroundfile

opts.trim = 64; %intial trim

opts.fineGrid = 64; %fine trim (min 256 for demonstrator b/s strehl computation)

opts.fitGrid = 64; %fine trim

opts.threshold = 10; % min counts threshold

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
filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\2017_07_08\DX_center_new\Spooled files.sifx';
% filename = 'P:\iLocater\NIC\2017_07_09\focus_25\Spooled files.sifx';
% filename = 'P:\iLocater\QuadCell\IR\timing\Spooled files.sifx';
% filename = 'P:\iLocater\QuadCell\IR\quad_cell\cx_04\Spooled files.sifx';
% filename = 'P:\iLocater\QuadCell\IR\Andor\sine_uncorrected\Spooled files.sifx';
% filename = 'C:\Data\2018_04_23\Quad Cell Response Time Testing 0.1V offset\quad_cell\sine_2\Spooled files.sifx';

% filename = 'C:\Data\2018_05_09\set1\Spooled files.sifx';
%filename = 'C:\Data\2018_05_14\set4_2\Spooled files.sifx';

% outfile = 'C:\Data\2018_05_09\set1\trace.dat';
test = Andor(filename,opts,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
toc
test = analyzeAndorData(test,opts); % reads in data and populates fundamental object properties
test.psfPlot
test.histPlot
toc
% save SXFast2 SXFast2
% 
% figure(10000);
% hold on
% trace(:,1) = test.time(:,1)*1e-6;
% trace(:,2) = test.fitParameters(:,2);
% % save(outfile,'trace','-ascii');
% plot(trace(:,1), trace(:,2));

% figure(2); 
% hold on; 
% plot(cx9.time(:,1)*1e-6,cx9.fitParameters(:,2))
% plot(cx8.time(:,1)*1e-6,cx8.fitParameters(:,2))
% plot(cx7.time(:,1)*1e-6,cx7.fitParameters(:,2))
% plot((cx04.time(:,1)*1e-6)-0.015,cx04.fitParameters(:,2))
% xlabel('time (s)')
% ylabel('x posistion (pix)')
% legend ('cx=9','cx=8','cx=7','cx=0.4')