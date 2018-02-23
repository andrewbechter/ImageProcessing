classdef Basler <Image
    properties
    end
    
    methods
        function [obj] = Basler(directory,memoryStep,startFrame,endFrame,type) % constructor (fills in field values with default settings)   
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                memoryStep = 10;
                type = 'fast';
                endFrame = 0;
                startFrame = 1;
            elseif nargin < 3
                type = 'fast';
                endFrame = 0;
                startFrame = 1;
            elseif nargin < 4
                type = 'fast';
                endFrame = 0;
            elseif nargin < 5
                type = 'fast';
            end
            
            delta = 40; % sets the frame size for fitting (2*delta)
            obj.memoryStep = memoryStep;
            obj.timeUnits = 1; %recorded in seconds
            obj.filename = directory;% stores the filename with the object
            
            %First get data needed on entire dataSet before reading in each
            %frame individually
            [obj.noFrames,allnames] = Basler.getDataSetInfo(directory); %check set exists, frame dimensions, total number of frames)
            
            jj = 0;% counter for tracking stored frames
            
            if endFrame > obj.noFrames || endFrame == 0
                endFrame = obj.noFrames;
            end
            
            for ii = startFrame:endFrame
                %read each frame individually
                %=============================================%
                %purpose: grabs a single frame and time stamp
                %inputs: directory for Basler image files
                %outputs: frame, time stamps
                [imageData] = Basler.getFrame(directory,allnames(ii));
                obj.dimensions = size(imageData);  % assign dimensions to object (used for fitting and plotting)
                %=============================================%
                %purpose: fastpeak finder and initial guesses for trim and fit
                %inputs: frame
                %outputs: initial values format[amp,cenX,sigmaX,cenY,sigmaY,theta,offset];
                [obj.iVals(ii,:),obj.totalCounts(ii),obj.flag(ii) ] = Image.xcorrFit(imageData);
                %[iVals(ii,:)] = Andor.coarsefit(imageData);
                
                if strcmp(type,'full')==1 && obj.flag(ii) == 0
                    
                    %=============================================%
                    %purpose: check boundaries of cuts vs. the edges of frame
                    %inputs: ycenter, xcenter, psf sigma, dimensions, delta
                    %outputs: cut locations for frame
                    [cuts(ii,:)] = Image.findFrameCuts(obj.iVals(ii,4),obj.iVals(ii,2),obj.dimensions,delta);
                    
                    
                    %=============================================%
                    %inputs: ycenter, xcenter, constant sigma, dimensions, delta
                    %outputs: cut locations for frame
                    [cutFrame,xdata,flag] = Image.trimFrame(imageData,cuts(ii,:));%trim the frame
                    
                    if flag ==1
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
                        fprintf('%i\n',ii)
                        obj.frame(:,:,jj) = imageData;
                        obj.storedNums(jj) = ii;
                        obj.cuts(jj,:) = cuts(ii,:);
                    end
                    
                end
            end
        end
        %=============================================%
        % methods that bundle analysis methods %
        function [obj] = analyzeBaslerData(obj)
            [obj] = calcMean(obj);
            [obj] = calcRMS(obj);
            [obj] = calcRange(obj);
%             [obj] = calcFrameRate(obj);
%             [obj] = calcFFT(obj);
%             [obj] = calcPSD(obj);
        end
        function [obj] = quickAnalyzeAndorData(obj)
        end
        %=============================================%
        % methods that use individual frames 
        function [obj] = inspectFrame(obj,frameNumber)
            [obj.noFrames,allnames] = Basler.getDataSetInfo(obj.filename); %check set exists, frame dimensions, total number of frames)
            [imageData] = Basler.getFrame(obj.filename,allnames(frameNumber));
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
                type = '.tiff'; % need to change this to a check
                fileID = strcat(directory,'/*',type);
                allFiles = dir(fileID);
                allNames = {allFiles.name};
                noFrames = size(allNames,2); %M is the number of files in the directory
        end
        function [imageData] = getFrame(path,name)
            filename = strcat(path,'/',name);
            imageData = double(imread(char(filename)));
        end
        %=============================================%
    end
end
