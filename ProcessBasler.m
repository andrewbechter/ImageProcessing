% ProcessBasler

memStep = 1; 
startFile = 1;
endFile = 0;
fs = 1;


filename = '/Volumes/Projects/iLocater/AcquisitionCamera/COTS/BackIllumination/Set1';
tic
set1 = Basler(filename,memStep,startFile,endFile,'full',fs); % reads in data and populates fundamental object properties
toc

set1 = analyzeBaslerData(set1); %calculates standard analysis
set1.psfPlot
set1.histPlot
set1.Mean(2)
set1.Mean(4)
