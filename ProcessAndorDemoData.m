%MainProcessingScript
clear
tic
%function descriptions:
%name: Andor
%purpose: grabs a single frame and time stamp
%inputs: ## filaname: filanem including directory and .sifx
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

memStep = 1000;
startFile = 1;
endFile = 0;

%=============================================%
%AUSTRALIS
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx';
% Aus1 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus1 = analyzeAndorData(Aus1); %calculates standard analysis
% save Aus1 Aus1
% clear Aus1 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_2\Spooled files.sifx';
% Aus2 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus2 = analyzeAndorData(Aus2); %calculates standard analysis
% save Aus2 Aus2
% clear Aus2 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_3\Spooled files.sifx';
% Aus3 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus3 = analyzeAndorData(Aus3); %calculates standard analysis
% save Aus3 Aus3
% clear Aus3 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_4\Spooled files.sifx';
% Aus4 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus4 = analyzeAndorData(Aus4); %calculates standard analysis
% save Aus4 Aus4
% clear Aus4 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_153_sub\Spooled files.sifx';
% Aus153 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus153 = analyzeAndorData(Aus153); %calculates standard analysis
% save Aus153 Aus153
% clear Aus153 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_5\Spooled files.sifx';
% Aus5 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus5 = analyzeAndorData(Aus5); %calculates standard analysis
% save Aus5 Aus5
% clear Aus5 filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_6\Spooled files.sifx';
% Aus6 = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Aus6 = analyzeAndorData(Aus6); %calculates standard analysis
% save Aus6 Aus6
% clear Aus6 filename

%=============================================%
%72LEO
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\72Leo\72Leo_full\Spooled files.sifx';
% Leo72Full = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Leo72Full = analyzeAndorData(Leo72Full); %calculates standard analysis
% save Leo72Full Leo72Full
% clear Leo72Full filename

% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\72Leo\72Leo_sub\Spooled files.sifx';
% Leo72Sub = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% Leo72Sub = analyzeAndorData(Leo72Sub); %calculates standard analysis
% save Leo72Sub Leo72Sub
% clear Leo72Sub filename

%=============================================%
%KXVIR
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\KXVIR\KXVIR_full\Spooled files.sifx';
% KXVIRFull = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% KXVIRFull = analyzeAndorData(KXVIRFull); %calculates standard analysis
% save KXVIRFull KXVIRFull
% clear KXVIRFull filename
% 
% filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\KXVIR\KXVIR_sub\Spooled files.sifx';
% KXVIRSub = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties
% KXVIRSub = analyzeAndorData(KXVIRSub); %calculates standard analysis
% save KXVIRSub KXVIRSub
% clear KXVIRSub filename
