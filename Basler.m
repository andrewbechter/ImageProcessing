classdef Basler <Image
    properties
        timeZero        % start of data set (may not need to keep this)
        focal_length    % focal length of the
        Fnum            % f number
        pix             % pixel size (meters)
        gridX           % processed frame grid
        gridY           % processed frame grid
        instSR          % instantaneous Strehl Ratio
        imStack         % stacked images (non centered)
        r               % 0-centered radius calculated from cenx, ceny
        allnames        % all file extension names in directory specified
    end 
    
    methods
        %--------------%
        % constructor
        %--------------%
        
        
        function [obj] = Basler(directory,opts) % constructor (fills in field values with default settings)
            
            %%----------------------------
            % Break out processing options
            %%----------------------------
            
            memoryStep = opts.memStep ;
            
            endFrame = opts.endFile;
            
            startFrame = opts.startFile;
            
            fs = opts.fs;
            
            type = opts.type ;
            
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                disp('you must specify processing options')              
            end
            
            %%----------------------------
            % Break out processing options
            %%----------------------------
            
            backgroundFile = opts.backgroundfile;
            
            trim = opts.trim;
            
            srgrid = opts.srgrid;
            
            fineGrid = opts.fineGrid;
            
            fitGrid = opts.fitGrid;
            
            obj.memoryStep = memoryStep;
            
            obj.timeUnits = 1; %recorded in seconds
            
            obj.filename = directory;% stores the filename with the object
            
            %%----------------------------
            % Break out optical parameters
            %%----------------------------
            
            lambda = opts.lambda; % wavelegnth in meters
            
            focal_length = opts.focal_length; % effective focal length
            
            d = opts.d; % beam diameter
            
            pix = opts.pix; %pixel size
            
            npup = opts.npup; %number of samples across the pupil
            
            alpha = opts.alpha; %blocking fraction
            
            F = focal_length/d; % compute F/number
            
            q = F*lambda/(pix); %samples across PSF
            
            %%----------------------------
            % Compute diffraction pattern
            %%----------------------------
            
            [airy] = Image.diffraction_pattern(q,lambda,npup,alpha); % compute diffraction pattern
            
            [airy,~,~] = Image.centroid_center(airy, 0.1, 5, false); % center image and calculate dx,dy
            
            obj. airy = airy;
            
            %%--------------------
            % Initialize variables
            %%--------------------
            
            obj.imStack =[];
            
            obj.coAdd=[];
            
            %%------------------
            % Check data set
            %%------------------
            %purpose:First get data needed on entire dataSet before reading in each frame individually
            %inputs:
            %outputs:
            
            [obj.noFrames,allnames] = Basler.getDataSetInfo(directory); %check set exists, frame dimensions, total number of frames)
            obj.allnames = allnames;
            obj.time = 0:1/fs:(obj.noFrames-1)/fs;% time vector is purely calculated based on sampling frequency (frame rate)
            
            obj.frameRate = fs; % attach frame rate in Hz to object.       

            
            %%--------------------
            % Assign set start/end
            %%--------------------
            
            jj = 0;% counter for tracking stored frames
            
            if endFrame > obj.noFrames || endFrame == 0
                endFrame = obj.noFrames;
            end
            
            chk1 = round(0.1*(endFrame-startFrame));
            chk2 = round(0.5*(endFrame-startFrame));
            chk3 = round(0.9*(endFrame-startFrame));
            
            
            %%-----------------------
            % Start processing frames
            %%-----------------------
            
            for ii = startFrame:endFrame
                
                if ii == startFrame+chk1
                    fprintf('Working on frame %i ...10%% done\n',ii)
                elseif ii == startFrame+chk2
                    fprintf('Working on frame %i ...50%% done\n',ii)
                elseif ii == startFrame+chk3
                    fprintf('Working on frame %i ...90%% done\n',ii)
                end
                                
                            
                %%-----------------------------
                % Read each frame (indiviually)
                %%-----------------------------
                %=============================================%
                %purpose: grabs a single frame and time stamp
                %inputs: directory for Basler image files
                %outputs: frame, time stamps
                [imageData] = Basler.getFrame(directory,allnames(ii));
                
                obj.dimensions = size(imageData);  % assign dimensions to object (used for fitting and plotting)
                
                height = obj.dimensions(1);  % assign dimensions to object (used for fitting and plotting)
                
                width = obj.dimensions(2);
                
                %%-----------------
                % Remove background
                %%-----------------
                
                [imageData,~,~] = Image.subBackground(imageData,backgroundFile);
                               
                %imageData(imageData<5) = 0;

                %%----------------
                % Constuctt Grid
                %%----------------
                
                sizex = size(imageData,2);
                
                sizey = size(imageData,1);
                
                [X,Y] = meshgrid(1:sizex,1:sizey);
                
                trim = min([trim,width,height]);
                
                X= Image.pad(X,trim,trim);
                
                Y= Image.pad(Y,trim,trim);
                
                %%------------------------
                % Center frames and trim
                %%------------------------
                
                [frame] = Image.pad(imageData,trim,trim);% standard size if frame is larger need offsets in x,y

                if isempty(obj.imStack) == 1
                    obj.imStack = frame; % stack and add frames
                else
                    obj.imStack = obj.imStack + frame; % stack and add frames
                end
                
                [frame,dx,dy] = Image.centroid_center(frame, 0.1, 5, false); % center image and calculate dx,dy
                
                X = X + dx; %adjust grid to be centered
                
                Y = Y + dy; %adjust grid to be centered
                
                fineGrid = min([fineGrid,width,height]);
                
                [frame] = Image.pad(frame,fineGrid,fineGrid);% standard size if frame is larger need offsets in x,y
                
                X= Image.pad(X,fineGrid,fineGrid); %trim grid
                
                Y= Image.pad(Y,fineGrid,fineGrid); %trim grid
                
                if max(max(frame)) < opts.threshold
                    obj.flag = 1;
                    
                    obj.fitParameters(ii,:) = [0,0,0,0,0,0,0];
                    
                    obj.instSR(ii) = 0;
                    
                else
                    
                    
                    if opts.fitPSF
                        %%------------------------
                        % Gaussian Fit
                        %%------------------------
                        
                        fitGrid = min([fitGrid,width,height]);
                        
                        [fitframe] = Image.pad(frame,fitGrid,fitGrid);% standard size if frame is larger need offsets in x,y
                        
                        xdata(:,:,1)= Image.pad(X,fitGrid,fitGrid); %trim grid
                        
                        xdata(:,:,2)= Image.pad(Y,fitGrid,fitGrid); %trim grid
                        
                        xcen = X(1,ceil(end/2));
                        
                        ycen = Y(ceil(end/2), 1);
                        
                        sig = size(fitframe)./5;
                        
                        amp = max(max(fitframe));
                        
                        obj.iVals(ii,:) = [amp,xcen,sig(1),ycen,sig(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
                        
                        [obj.fitParameters(ii,:)] = Image.subPixelPeakV3(fitframe,obj.iVals(ii,:),0,xdata);% sub pixel fitting
                    end
                    
                    %%--------------------
                    % Co-Add frames
                    %%--------------------
                    
                    if isempty(obj.coAdd) == 1
                        obj.coAdd = frame; % stack and add frames
                    else
                        obj.coAdd = obj.coAdd + frame; % stack and add frames
                    end
                    
                    if opts.instSR
                        %%--------------------
                        % Compute Strehl Ratio
                        %%--------------------
                        sr_frame = Image.pad(frame,srgrid,srgrid);
                        obj.instSR(ii) = Image.calculate_SR(sr_frame,airy,false);
                        
                    end
                end
                %%-------------
                % Saved Data
                %%-------------
                
                if(mod(ii,obj.memoryStep) == 0) || ii == 1 % store every nth frame and the 1st one
                    jj = jj+1;
                    obj.frame(:,:,jj) = frame;
                    obj.storedNums(jj) = ii;
                    obj.gridX(:,:,jj) = X;
                    obj.gridY(:,:,jj) = Y;
                end
            end            
        end
        
        
        function [obj] = old_Basler(directory,memoryStep,startFrame,endFrame,type,fs) % constructor (fills in field values with default settings)
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                memoryStep = 10;
                type = 'fast';
                endFrame = 0;
                startFrame = 1;
                fs = 1;
            elseif nargin < 3
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
            
            delta = 50; % sets the frame size for fitting (2*delta)
            obj.memoryStep = memoryStep;
            obj.timeUnits = 1; %recorded in seconds
            obj.filename = directory;% stores the filename with the object
            
            
            %First get data needed on entire dataSet before reading in each
            %frame individually
            [obj.noFrames,allnames] = Basler.getDataSetInfo(directory); %check set exists, frame dimensions, total number of frames)
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
                %inputs: directory for Basler image files
                %outputs: frame, time stamps
                [imageData] = Basler.getFrame(directory,allnames(ii));
                obj.dimensions = size(imageData);  % assign dimensions to object (used for fitting and plotting)
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
        function [obj] = analyzeBaslerData(obj,opts)
            [obj] = calcMean(obj);
            [obj] = calcRMS(obj);
            [obj] = calcRange(obj);
            obj.r = sqrt((obj.fitParameters(:,2)-obj.Mean(2)).^2 + (obj.fitParameters(:,4)-obj.Mean(4)).^2);
            [obj] = calcStrehl(obj,opts);
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
            imageData = double(imread(char(filename),'tiff'));
            bit_scale = min(diff(unique(imageData)));
            imageData = imageData./bit_scale;
        end
        %=============================================%
    end
end
