# ImageProcessing
Reads in frames from ANDOR Zyla and Basler detectors. Calculates centroids, 2D Gaussian fit parameters, FFTs, PSDs, and scatter statistics

### Prerequisites
Matlab and Matlab curve fitting toolbox is required. 
For native ANDOR files (.sifx, .sif), the MATLAB SIF reader is required.

## Processing a data set

The default class constructor methods "Basler" and "Andor" require the directory of the data set for processing. The default mode for data processing will start at the first file/frame and end with the last file/frame. Only the first frame is stored in the object variable to conserve memory. A full 2D Gaussian fit is performed on every frame. 

For BASLER files the filename syntax is as follows:
```
filename = 'filename = 'P:\iLocater\QuadCell\IR\Set118\'
labSet = Basler(filename);
```
For ANDOR files the filename syntax is as follows:
```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename);
```
The Australis variable will be a type Andor and the labSet variable will be type Basler. 

There are four options availible for processing data in different ways for Andor frames: stored frames, start frame , end Frame and fit type. For example:

```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename,storedFrames,startFrame,endFrame,fitType);
```

user inputs:
- 'storedFrame' can be any integer, n. This variable will store every nth frame, incuding the first. 0 will only store the first frame.
- 'startFrame' allows the user to adjust which frame to start in the set. 
- 'endFrame' allows the user to adjust the final frame. a zero value will process to the end. 
- 'fitType' can be either 'full' or 'fast'. 'full' uses the 2D Gausssian fit whereas 'fast' computes 1D slice fits about the cener pixel. Contrary to the name, 'fast' does not actually work that much faster so 'full' is recommended. Additional options will be present in the future for that. 

An additional input is required for Basler data, called sampling Frequency. This is used to calculate the time each frame was taken and defaults to 1Hz if not specified. These are only approximate values and will not be correct if the detector drops frames while recording!

Note the filename only needs to contain the directory (no file or extension) for Basler data. 

```
filename = 'P:\iLocater\QuadCell\BaslerSet1\'
Set1 = Basler(filename,storedFrames,startFrame,endFrame,fitType,samplingFrequency);
```

An example of a custom process run with Andor data:

```
Australis = Andor(filename,100,5,0,'full'); 
% Australis variable will store every 100th frame, starting with thr 5th and ending with the final frame using the 2D fit algorithm
```
The object files created after processing have the following object parameters. These parameters behave like standard variables and can be accessed using the dot notation. 

```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename,100,5,0,'full'); 
total frames = Austiralis.noFrames  % This command will retrieve the number of frames from the parameter noFrames (I know not a good name) and assign the value to total frames.
```
The dot notation is useful for checking parameters quickly in the command line or in an analysis script. 

## Analysing a data set

Most, if not all are parameters are populated after running the Andor/Basler function followed by the analysis function. 

```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename,100,5,0,'full'); 
Australis = analyzeAndorData(Austalis); % analysis
```
The analysis function runs through a preset list of analysis calculations using default settings. The only input required is the object (i.e. 'Australis'). The output must be specified to update or overwrite parameters in the object' Austalis'. If no output is assigned the result of the funtion will be created in 'ans' just like any other function in Matlab with no specified output variable. 

## Saving a data set

Saving can be done just like any other variable in matlab
```
save Australis Australis % saves the variable Australis with the .mat name, 'Australis'
```

## Loading a data set

Loading can be done just like any other variable in matlab
```
load('Australis') 
```
The path to the Class files must be in matlabs working directory to load a file correctly. This can be done by running the script add path, or by manually adding '/Volumes/Software/ImageProcessing' for mac and 'S:\ImageProcessing' for windows. 


The list of parameters in Andor/Basler is below:

```
        filename        % filename for later reference
        timeUnits       % e.g. Andor is recorderd in microseconds (1e-6 = seconds)
        frame           % stored frames
        storedNums      % stored frame numbers
        noFrames        % total frames in data set
        dimensions      % width x height
        time            % relative time (microseconds for Andor), absolute time (matlab), flag (0 = measured from stamps, 1 = calculated from kinetic cycle time)
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

```

The default analysis can be overridden by using individual functions:

Full method list for image processing class. These functions all require objects as inputs and outputs. Some optional inputs 
```
%=============================================%
% methods that calculate things using time series data from object %
function [obj] = calcFrameRate(obj) % reads in a frame, calculates frame rate from time stamp vecotor 
function [obj] = calcFFT(obj,xdata) % calculates single sided fourier transform. default is y position, xdata will override this 
function [obj] = calcPSD(obj,xdata) % calculates the power spectral density. default is y position, xdata will override this

%=============================================%
% methods that calculate image quality data from object %
function [obj] = calcCircPSF(obj,sigma) %idenitfy psfs with circular core, with a xsigma, ysigma smaller than input sigma. 
function [obj] = calcWFE(obj) % does nothing right now.
function [obj] = calcStrehlRatio(obj) % does nothing right now
function [obj] = calcFiberCoupling(obj) % does nothing right now

%=============================================%
% methods that calculate statistics data from fitParameters. These will exclude frames with 0 amplitude or bad frame flag 
function [obj] = calcMean(obj)
function [obj] = calcRange (obj)
function [obj] = calcRMS (obj)
function [obj] = calcDelta(obj) % delta is the radial distance from the mean position in x,y. 

%=============================================%
% methods that plot data from object
function psfPlot(obj,data_number)% default data_number is 1. User can specify any of the stored frames instead.
function histPlot(obj,x,y,stats)% default is to create x and y scatter historgams and stats. User can specify alternative data set
function FFTPlot(obj) % plots the FFT property in the object
function PSDPlot(obj) % plots the PSD property in the object
function inspectFrame(index) % plot any frame from raw data using the frame index

%=============================================%
% methods to save object%
function saveToStruct(obj, filename) % save the object as a default matlab structure. This remove the class/object nature from varaible
```

## Authors
* **Andrew Bechter** 

+See also the list of [contributors](https://github.com/andrewbechter/ImageProcessing/contributors) who participated in this project.

## License

...

## Acknowledgments

* Hat tip to anyone who's code was used...
* Fourier filter code directly used from T. C. O'Haver (toh@umd.edu),  version 1.5, May, 2007 
* Plotting inspired by Gero Nootz 2010 Fit 2D gaussian function to data
