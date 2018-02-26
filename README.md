# ImageProcessing
Reads in frames from ANDOR Zyla and Basler detectors. Calculates centroids, 2D Gaussian fit parameters, FFTs, PSDs, and scatter statistics

Processing ANDOR data ('.sifx' or '.sif') can be done using the command line or the 'Process.m' template script.
An example command: 

filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx';
Australis = Andor(filename); % reads in data and populates fundamental object properties using default methods and steps

function name: Andor ;
purpose: grabs a frame or frames from the .sifx spool, performance a basic fit, a full 2D fit, and retrieves the time stamp. 
inputs: filaname: filename including directory and .sifx (This input is required)
        memoryStep: interval to store frames (Total frames stored should be ~100 max, storing too many (10K) will use more memory). 
        frameStart: frame number to start on        
        frameEnd: frame number to end on 
        type: 'fast' or 'full' processing. 'fast' only returns approximate fit parameters, 'full' uses 2D Gaussian fit
        
notes: only the filename needs to be specified. The defaults for the rest of the inputs are : no frames
stored, start frame is 1, end frame is the last frame, and all frames processed using 'full'
outputs: Andor object with fit parameters, time stamps and stored frames if specified. 

memStep = 1000; %store every 1000th frame except for the first interval (Starts at 1, then 1000, 2000 etc.)
startFile = 1; % start at the first frame
endFile = 0; % go until the end 

filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx';
Australis = Andor(filename,memStep,startFile,endFile,'full'); % reads in data and populates fundamental object properties

Loading ANDOR variables:

All ANDOR data has been reduced and saved using the ImageProcessing software on the shared storage. The matlab data types are objects created by custom classes found in Software:\ImageProcessing\Andor. To load the files in matlab correctly, the above path to the Andor folder must be in your matlab path! If the file loads correclty the variable will immediately appread in the workspace as a type 'Andor'.

Once the files are loaded and recognized by matlab, object properties/methods contained in class files can be accessed using the dot notation.

For example:
"Australis.mat" is the stored matlab variable name.

Add "Software:\ImageProcessing\Andor" to path. You can do this by manually navigating and adding to path or running script "addPath.m"

"load(Australis.mat)" in the command line or double click the variable in an explorer window.
The object can be seen in the workspace. Double clicking the variable will reveal information/parameters avaible to the user. The same can be done by typing "Australis' into the terminal.

Analysing ANDOR variables:

%function name: analyzeAndorData
%purpose: calculates useful things using time stamps and fit parameters
%inputs: ## Andor object
%notes: this method modifies exisiting Andor objects by populating or
%overwriting useful things for analysis 
%outputs: ## Andor object


Parameters can be accessed just like standard variables when using the dot notation. 

For example:

"Australis.fitParameters" will show all the 2D fit values

Some built in functions can do statistical analysis and plotting
 
"Australis.psfPlot" will plot the first frame and show the fit values 

The property list for Andor Class is provided below 

        filename        % filename for later reference
        timeUnits       % Andor is recorderd in microseconds (time*1e-6 = seconds)
        frame           % stored frames
        storedNums      % stored frame numbers
        noFrames        % total frames in data set
        dimensions      % width x height
        time            % relative time(microseconds), absolute time (matlab), flag (0 = measured from stamps, 1 = calculated from kinetic cycle time)
        memoryStep      % frame storage interval
        iVals           % starting points for fit parameters. format x0 = [Amp,xo,wx,yo,wy,theta,offset];
        fitParameters   % calculated fit parameters. format x = [Amp,xo,wx,yo,wy,theta,offset];
        cuts            % locations to trim the frame
        frameRate       % calculated frame rate based off realtive time stamps (time(:,1))
        FFT             % calculates the single sided fast fourier transform based on sampling frequency
        PSD             % calculates the PSD based on sampling frequency
        Mean            % calculates the mean value of all fit parameters
        RMS             % calculates the RMS value of all fit parameters
        Range           % calculates the Range (max) value of all fit parameters
        totalCounts     % total number of counts in the whole frame
        flag            % if flag =1, frame trimming/ centroid failed, no 2D fit recorded.(i.e. horrible readout noise)
        circPSF         % index value of all the PSF's identified as small and circular (i.e. have a core)
        delta           % delta is the radial distance from the mean centroid.
        exFrames        % index values for frames which meet a certain condition and have been saved
        tempFrame       % temporary frames read in by the user for inspection

For a full list of functions see the better readMe page (working on it right now). 
