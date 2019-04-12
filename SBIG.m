classdef SBIG <Image
    properties
    end
    
    methods
        function [obj] = SBIG(directory,memoryStep,startFrame,endFrame,type,fs) % constructor (fills in field values with default settings)
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                memoryStep = 10;
                type = 'fast';
                endFrame = 0;
                startFrame = 1;
                fs = 1;
            elseif nargin < 1
                type = 'fast';
                endFrame = 0;
                startFrame = 1;
                fs = 1;
            elseif nargin < 4
                type = 'fast';
                endFrame = 0;
                fs= 1;
            elseif nargin < 5
                type = 'fast';
                fs = 1;
            elseif nargin < 6
                fs = 1;
                
            end
            detType = 'SBIG';
            delta = 50; % sets the frame size for fitting (2*delta)
            obj.memoryStep = memoryStep;
            obj.timeUnits = 1; %recorded in seconds
            obj.filename = directory;% stores the filename with the object
            
            
            %First get data needed on entire dataSet before reading in each
            %frame individually
            [obj.noFrames,allnames] = SBIG.getDataSetInfo(directory); %check set exists, frame dimensions, total number of frames)
            obj.time = 0:1/fs:(obj.noFrames-1)/fs;% time vector is purely calculated based on sampling frequency (frame rate)
            obj.frameRate = fs; % attach frame rate in Hz to object.
            jj = 0;% counter for tracking stored frames
            
            if endFrame > obj.noFrames || endFrame == 0
                endFrame = obj.noFrames;
            end
            
            chk1 = round(0.1*(endFrame-startFrame));
            chk2 = round(0.5*(endFrame-startFrame));
            chk3 = round(0.9*(endFrame-startFrame));
            
            for ii = startFrame:endFrame
                obj.flag(ii) =0;
                %read each frame individually
                if ii == startFrame+chk1
                    fprintf('Working on frame %i ...10%% done\n',ii)
                elseif ii == startFrame+chk2
                    fprintf('Working on frame %i ...50%% done\n',ii)
                elseif ii == startFrame+chk3
                    fprintf('Working on frame %i ...90%% done\n',ii)
                end
                
                %=============================================%
                %purpose: grabs a single frame and time stamp
                %inputs: directory for SBIG image files
                %outputs: frame, time stamps
                [imageData] = SBIG.getFrame(directory,allnames(ii));
                
                obj.dimensions = size(imageData);  % assign dimensions to object (used for fitting and plotting)
                
                background = fitsread('S:\ImageProcessing\CalibrationFiles\SBIG\Background\COTS.fit');
                [imageData] = Image.subBackground(imageData,background./1.25); %remove brackground from SBIG detector
                [imageData] = Image.sigmaClipping(imageData,background);
                %=============================================%
                %purpose: fastpeak finder and initial guesses for trim and fit
                %inputs: frame
                %outputs: initial values format[amp,cenX,sigmaX,cenY,sigmaY,theta,offset];
                [obj.iVals(ii,:),obj.totalCounts(ii),obj.flag(ii) ] = Image.xcorrFit(imageData);
                %
                if ii == startFrame
                    
                    [obj.flag(ii),sigma]= Image.fit1D(obj.iVals(ii,1),[obj.iVals(ii,2),obj.iVals(ii,4)],imageData);
                    %[iVals(ii,:)] = Andor.coarsefit(imageData);
                end
                
                obj.iVals(ii,3) = sigma(1); % use startFrame sigma as guess for all sigmas
                obj.iVals(ii,5) = sigma(2); % use startFrame sigma as guess for all sigmas
                
                %=============================================%
                %purpose: check boundaries of cuts vs. the edges of frame
                %inputs: ycenter, xcenter, psf sigma, dimensions, delta
                %outputs: cut locations for frame
                [cuts(ii,:)] = Image.findFrameCuts(obj.iVals(ii,4),obj.iVals(ii,2),obj.dimensions,delta);
                
                if strcmp(type,'full')==1
                    
                    %=============================================%
                    %inputs: ycenter, xcenter, constant sigma, dimensions, delta
                    %outputs: cut locations for frame
                    [cutFrame,xdata,flag] = Image.trimFrame(imageData,cuts(ii,:));%trim the frame
                    
                    if obj.flag ==1
                        obj.fitParameters(ii,:) = zeros(1,length(obj.iVals(ii,:)));
                        obj.flag(ii) = 1;
                    else
                        %=============================================%
                        %inputs: ycenter, xcenter, constant sigma, dimensions, delta
                        %outputs: cut locations for frame
                        [obj.fitParameters(ii,:)] = Image.subPixelPeakV3(cutFrame,obj.iVals(ii,:),0,xdata);% sub pixel fitting
                    end
                    
                    
                    if(mod(ii,obj.memoryStep) == 0) || ii == 1 % store every nth frame and the 1st one
                        jj = jj+1;
                        %fprintf('%i\n',ii)
                        obj.frame(:,:,jj) = imageData;
                        obj.storedNums(jj) = ii;
                        obj.cuts(jj,:) = cuts(ii,:);
                    end
                end
            end
        end
        %=============================================%
        % methods that bundle analysis methods %
        function [obj] = analyzeSBIGData(obj)
            [obj] = calcMean(obj);
            [obj] = calcRMS(obj);
            [obj] = calcRange(obj);
            [obj] = addFrames(obj);
            %             [obj] = calcFrameRate(obj);
            %             [obj] = calcFFT(obj);
            %             [obj] = calcPSD(obj);
        end
        function [obj] = quickAnalyzeSBIGData(obj)
        end
        %=============================================%
        % methods that use individual frames
        function [obj] = inspectFrame(obj,frameNumber)
            [obj.noFrames,allnames] = SBIG.getDataSetInfo(obj.filename); %check set exists, frame dimensions, total number of frames)
            [imageData] = SBIG.getFrame(obj.filename,allnames(frameNumber));
            obj.tempFrame = imageData;
            figure
            imagesc(obj.tempFrame)
        end
        %=============================================%
    end
    
    methods(Static)
        %=============================================%
        % functions that recover the image and time stamps %
        function [noFrames,allNames] = getDataSetInfo(directory)
            type = '.fit'; % need to change this to a check
            fileID = strcat(directory,'/*',type);
            allFiles = dir(fileID);
            allNames = {allFiles.name};
            noFrames = size(allNames,2); %M is the number of files in the directory
        end
        function [imageData] = getFrame(path,name)
            filename = strcat(path,'/',name);
            imageData = double(fitsread(char(filename)));
        end
        %=============================================%
    end
end
