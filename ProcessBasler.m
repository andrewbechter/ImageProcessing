% ProcessBasler

memStep = 100; 
startFile = 1;
endFile = 0;
fs = 1;

% filename = '/Volumes/Projects/iLocater/AcquisitionCamera/COTS/FiberChannel/Set17/';
% filename = 'P:\iLocater\QuadCell\IR\Set160\';
filename = 'P:\iLocater\WhitelaseCoupler\Set7\';
tic
set1 = Basler(filename,memStep,startFile,endFile,'full',fs); % reads in data and populates fundamental object properties
toc

set1 = analyzeBaslerData(set1); %calculates standard analysis
set1.psfPlot
set1.histPlot
set1.Mean(2)
set1.Mean(4)
