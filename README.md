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

There are four options availible for processing data in different ways: stored frames, start frame , end Frame and fit type. 

```
filename = 'P:\iLocater\iLocater_Demonstrator\LBT_Data\Forerunner\2016_04_18\Australis\Australis_1\Spooled files.sifx'
Australis = Andor(filename,storedFrames,startFrame,endFrame,fitType);
```

storedFrame can be any integer, n. This variable will store every nth frame, incuding the first. 0 will not store any frames
startFrame allows the user to adjust which frame to start in the set. 


## Authors

* **Andrew Bechter** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

N/A

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc




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
