# ImageProcessing
Reads in frames from ANDOR Zyla and Basler detectors. Calculates centroids, 2D Gaussian fit parameters, FFTs, PSDs, and scatter statistics

Processing ANDOR data



Loading ANDOR variables

All ANDOR data has been reduced and saved using the ImageProcessing software on the shared storage. 
The matlab data types are objects created by custom classes found in Software:\ImageProcessing\Andor. 
To load the files in matlab correctly, the above path to the Andor folder must be in your matlab path!
Once the files are loaded and recognized by matlab, object properties/methods contained in class files can be accessed using the dot notation

For example:

"Australis.mat" is the stored matlab variable name

Add "Software:\ImageProcessing\Andor" to path (you can do this by manually navigating and adding to path or running script "addPath.m"

Type "load(Australis.mat)" in the command line or double click the variable in an explorer window

The object can be seen in the workspace. Double clicking the variable will reveal information/parameters avaible to the user. The same can be done by typing "Australis' into the terminal.

Parameters can be accessed just like standard variables when using the dot notation. 

For example:

"Australis.fitParameters" will show all the 2D fit values

Some built in functions can do statistical analysis and plotting
 
"Australis.psfPlot" will plot the first frame and show the fit values 


For a full list of functions see 
