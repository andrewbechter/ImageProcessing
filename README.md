# ImageProcessing
Reads in frames from ANDOR Zyla and Basler detectors. Calculates centroids, 2D Gaussian fit parameters, FFTs, PSDs, and scatter statistics

### Prerequisites
Matlab and Matlab curve fitting toolbox is required. 
For native ANDOR files (.sifx, .sif), the MATLAB SIF reader is required.

## Processing a data set

The default class constructor methods "Basler" and "Andor" require the directory of the data set for processing. The default mode for data processing will start at the first file/frame and end with the last file/frame. Only the first frame is stored in the object variable to conserve memory. A full 2D Gaussian fit is performed on every frame. 

For BASLER files the filename synatax is as follows:
```
filename = 'filename = 'P:\iLocater\QuadCell\IR\Set118\'
labSet = Basler(filename);
```
For ANDOR files the filename synatax is as follows:
```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename);
```
The Australis variable will be a type Andor and the labSet variable will be type Basler. 

There are four options availible for processing data in different ways: stored frames, start frame , end Frame and fit type. For example:

```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename,storedFrames,startFrame,endFrame,fitType);
```
user inputs:
- 'storedFrame' can be any integer, n. This variable will store every nth frame, incuding the first. 0 will only store the first frame.
- 'startFrame' allows the user to adjust which frame to start in the set. 
- 'endFrame' allows the user to adjust the final frame. a zero value will process to the end. 
- 'fitTpye' can be either 'full' or 'fast'. 'full' uses the 2D Gausssian fit whereas 'fast' computes 1D slice fits about the cener pixel. Contrary to the name, 'fast' does not actually work that much faster so 'full' is recommended. Additional options will be present in the future for that. 

```
Australis = Andor(filename,100,5,0,'full'); 
% Australis variable will store every 100th frame, starting with thr 5th and ending with the final frame using the 2D fit algorithm
```
object parameters
```
        filename        % filename for later reference
        timeUnits       % e.g. Andor is recorderd in microseconds (1e-6 = seconds)
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

```
For a full list of functions see a better README (working on it right now). 

## Authors
* **Andrew Bechter** 

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

N/A

## Acknowledgments

* Hat tip to anyone who's code was used (Original plotting code, Gaussian 2D, ifilter)
* Inspiration
* etc
