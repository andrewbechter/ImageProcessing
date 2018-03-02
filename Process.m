%Process Image Data
memStep = 1000; % 'first'
startFile = 1;
endFile = 0;

tic
filename = 'P:\iLocater\NIC\2017_07_09\SX_fast3\Spooled files.sifx';
SXFast3 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
toc
SXFast3 = analyzeAndorData(SXFast3); % reads in data and populates fundamental object properties
SXFast3.psfPlot
SXFast3.histPlot
toc














