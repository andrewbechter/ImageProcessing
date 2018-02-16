%MainProcessingScript
clear
tic
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';


Aus6 = Andor(filename,1000,[],'full'); % reads in data and populates fundamental object properties
toc
Aus6 = analyzeAndorData(Aus6); %calculates standard analysis
toc 
% psfPlot(Aus6)
% histPlot(Aus6)
% FFTPlot(Aus6)
% PSDPlot(Aus6)

