%Process Image Data
memStep = 10;
startFile = 1;
endFile = 0;

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx';
% set = Andor(filename,memStep,1,100,'full'); % reads in data and populates fundamental object properties

filename = 'P:\iLocater\QuadCell\IR\Set95\';
set = Basler(filename,memStep,1,100,'full'); % reads in data and populates fundamental object properties
set = analyzeBaslerData(set); %calculates standard analysis
set.psfPlot
set.histPlot

