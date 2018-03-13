%Process Image Data
memStep = 1000; % 'first'
startFile = 1;
endFile = 0;

tic
filename = 'P:\iLocater\NIC\2017_07_07\SX_fast2\Spooled files.sifx';
SXFast2 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
toc
SXFast2 = analyzeAndorData(SXFast2); % reads in data and populates fundamental object properties
SXFast2.psfPlot
SXFast2.histPlot
toc

save SXFast2 SXFast2