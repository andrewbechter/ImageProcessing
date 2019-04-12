classdef Andor < Image
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
        drift_x         % mean drift in x and y 
        drift_y 
        drift           % rss of x and y (pythagoras)
    end
    
    methods
        
        %--------------%
        % constructor
        %--------------%
        
        function [obj] = Andor(andor_file,opts)
            
            %%----------------------------
            % Break out processing options
            %%----------------------------
            
            memoryStep = opts.memStep ;
            
            endFrame = opts.endFile;
            
            startFrame = opts.startFile;
            
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
            
            obj.filename = andor_file;% stores the filename with the object
            
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
            
            [signal,obj.noFrames,and_size,width,height] = Andor.getDataSetInfo(andor_file); %check set exists, frame dimensions, total number of frames)
            
            obj.dimensions = [height,width];  % assign dimensions to object (used for fitting and plotting)
            
            [obj.timeZero] = Andor.getTimeZero (signal); % get the initial time (used to build computer clock time vector)
            
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
                
                %purpose: grabs a single frame and time stamp
                %inputs: flag, size of linearized frame, frame width, frame height, frame number
                %outputs: frame, time stamps
                
                [imageData,obj.time(ii,:)] = Andor.getFrameInfo(signal,and_size,width,height,ii);
                
                %%-----------------
                % Remove background
                %%-----------------
                
                [imageData,obj.medback,obj.sigback] = Image.subBackground(imageData,backgroundFile);
                
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
                
                X = X + dx; %adjust grid to account for PSF centering
                
                Y = Y + dy; %%adjust grid to account for PSF centering
                
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
            
            obj.time(:,2) = obj.time(:,2)+ obj.timeZero; % add initial time to matlab time to get the computer time
            
        end
        
        %------------------------------%
        % methods that bundle analysis
        %------------------------------%
        
        function [obj] = analyzeAndorData(obj,opts)
            
            [obj] = calcMean(obj);
            [obj] = calcRMS(obj);
            [obj] = calcRange(obj);
            [obj] = FrameRate(obj);
            obj.r = sqrt((obj.fitParameters(:,2)-obj.Mean(2)).^2 + (obj.fitParameters(:,4)-obj.Mean(4)).^2);
            [obj] = calcFFT(obj);
            [obj] = calcPSD(obj);
            [obj] = calcStrehl(obj,opts);
            %[obj] = calcCentroidDrift(obj);
            
        end
        
        function [obj] = quickAnalyzeAndorData(obj)
        end
        
        %------------------------------%
        % Andor Analysis
        %------------------------------%
        
        function [obj] = calcCentroidDrift(obj)
            %function descriptions:
            
            %name: calcCentroidDrift
            %purpose: linear fit to the x and y centroid 
            %inputs: ANDOR object
            %outputs: fit parameters measured in pixels         
            
            %%--------------------------
            % Fit X and Y centroids
            %%--------------------------
            
            x =  obj.fitParameters(:,2);
            y =  obj.fitParameters(:,4);
            ind1 = x>mean(x)-3*std(x(1:1000));
            ind2 = x<mean(x)+3*std(x(1:1000));
            indx = logical(ind1.*ind2);
            
            ind1 = y>mean(y)-3*std(y(1:1000));
            ind2 = y<mean(y)+3*std(y(1:1000));
            indy = logical(ind1.*ind2);
            
            ind = logical(indx.*indy);
            x = x(ind);
            y = y(ind);
            t =  obj.time(:,1)*1e-6; % convert time into seconds
            t=t(ind);
            
            [px,Sx] = polyfit(t,x,1);
            [py,Sy] = polyfit(t,y,1);
            
            %%--------------------------
            % Assign to object
            %%--------------------------
            
            obj.drift_x = px(1);
            obj.drift_y = py(1);
            obj.drift = sqrt(px(1).^2+py(1).^2);
            
        end
        %-----------------------%
        % plot data from object
        %-----------------------%
        function [obj] = inspectFrame(obj,frameNumber)
            [signal,~,and_size,width,height] = Andor.getDataSetInfo(obj.filename);
            [imageData,~] = Andor.getFrameInfo(signal,and_size,width,height,frameNumber);
            obj.tempFrame = imageData;
            figure
            imagesc(obj.tempFrame)
        end
        
    end
    methods(Static)
        %---------------------------------%
        %Recover the image and time stamps
        %---------------------------------%
        function [timeZero] = getTimeZero (signal)
            [rc,format_time] = atsif_getpropertyvalue(signal,'FormattedTime');
            frame_info{1,1} = 'Formatted Time';
            frame_info{1,2} = format_time;
            timeZero = datenum(format_time,'ddd mmm dd HH:MM:SS yyyy');
        end
        
        function [imageData,frame_time] = getFrameInfo(signal,and_size,width,height,frame_number)
            
            timeflag = 0;
            % Below grabs a Frame
            [rc,data]=atsif_getframe(signal,frame_number-1,and_size);
            [rc,pattern]=atsif_getpropertyvalue(signal,'ReadPattern');
            
            if(pattern == '0') % not sure what this means but it is not for standard images
                calibvals = zeros(1,and_size);
                for i=1:and_size
                    [rc,calibvals(i)]=atsif_getpixelcalibration(signal,xaxis,(i));
                end
                plot(calibvals,data);
                title('spectrum');
                [rc,xtype]=atsif_getpropertyvalue(signal,'XAxisType');
                [rc,xunit]=atsif_getpropertyvalue(signal,'XAxisUnit');
                [rc,ytype]=atsif_getpropertyvalue(signal,'YAxisType');
                [rc,yunit]=atsif_getpropertyvalue(signal,'YAxisUnit');
                xlabel({xtype;xunit});
                ylabel({ytype;yunit});
                
            elseif(pattern == '4') % this check seems to be for images
                
                %This shapes a frame
                imageData = double(reshape(data,width,height));
                imageData = rot90(imageData);
                
                prop_str = ['TimeStamp ',num2str(frame_number)]; %prep time stamp variable
                [rc,time] = atsif_getpropertyvalue(signal,prop_str); %grab the time stamp
                
                if(frame_number > 1 && str2num(time) == 0) % if there is no time stamp (i.e == 0)
                    prop_str = 'KineticCycleTime';%prep cycle time variable
                    [rc,time_raw] = atsif_getpropertyvalue(signal,prop_str); %grab cycle time variable
                    time = (str2num(time_raw) * (frame_number-1)) * 1e6; % calculate time based on cycle time (in microseconds)
                    frame_time(1,1) = time; % store calculated time in frame_time array
                    timeflag = 1;
                else
                    frame_time(1,1) = str2num(time); % store the measured time vector in frame_time array
                end
                frame_time(1,2) = frame_time(1,1) / (24 * 60 * 60 * 1e6) - (10/24);% calculate the matlab time
                frame_time(1,3) = timeflag;
            end
        end
        
        %-------------%
        %Data set info
        %-------------%
        
        function [signal,no_frames,and_size,width,height] = getDataSetInfo(andor_file)
            rc=atsif_setfileaccessmode(0);
            rc=atsif_readfromfile(andor_file);
            
            if (rc == 22002)
                signal=0;
                [rc,present]=atsif_isdatasourcepresent(signal);
                if present
                    [rc,no_frames]=atsif_getnumberframes(signal);
                    if (no_frames > 0)
                        [rc,and_size]=atsif_getframesize(signal);
                        [rc,left,bottom,right,top,hBin,vBin]=atsif_getsubimageinfo(signal,0);
                        xaxis=0;
                        width = ((right - left)+1)/hBin;
                        height = ((top-bottom)+1)/vBin;
                    end
                end
            end
        end
        
        %----------------------------------------------%
        %old function for reading andor data (not used)%
        %----------------------------------------------%
        
        function [image_data] = readData(andor_file)
            image_data = [];
            frame_limit = 0; %If 0 or less, reads all frames. Otherwise only reads this many frames
            rc=atsif_setfileaccessmode(0);
            rc=atsif_readfromfile(andor_file);
            
            if (rc == 22002)
                signal=0;
                [rc,present]=atsif_isdatasourcepresent(signal);
                if present
                    [rc,no_frames]=atsif_getnumberframes(signal);
                    if(frame_limit > 0)
                        no_frames = frame_limit;
                    end
                    if (no_frames > 0)
                        [rc,and_size]=atsif_getframesize(signal);
                        [rc,left,bottom,right,top,hBin,vBin]=atsif_getsubimageinfo(signal,0);
                        xaxis=0;
                        width = ((right - left)+1)/hBin;
                        height = ((top-bottom)+1)/vBin;
                        frame_info = cell(1,2);
                        
                        for ii = 1: no_frames
                            % Below grabs a Frame
                            [rc,data(:,ii)]=atsif_getframe(signal,ii-1,and_size);
                            [rc,pattern]=atsif_getpropertyvalue(signal,'ReadPattern');
                            
                            if(pattern == '0') % not sure what this means but it is not for standard images
                                calibvals = zeros(1,and_size);
                                for i=1:and_size,[rc,calibvals(i)]=atsif_getpixelcalibration(signal,xaxis,(i));
                                end
                                plot(calibvals,data);
                                title('spectrum');
                                [rc,xtype]=atsif_getpropertyvalue(signal,'XAxisType');
                                [rc,xunit]=atsif_getpropertyvalue(signal,'XAxisUnit');
                                [rc,ytype]=atsif_getpropertyvalue(signal,'YAxisType');
                                [rc,yunit]=atsif_getpropertyvalue(signal,'YAxisUnit');
                                xlabel({xtype;xunit});
                                ylabel({ytype;yunit});
                                
                            elseif(pattern == '4') % this check seems to be for images
                                
                                %                                 This shapes a frame
                                width = ((right - left)+1)/hBin;
                                height = ((top-bottom)+1)/vBin;
                                newdata = double(reshape(data,width,height));
                                newdata = flipud(rot90(newdata));
                                
                                prop_str = ['TimeStamp ',num2str(ii)]; %prep time stamp variable
                                [rc,time] = atsif_getpropertyvalue(signal,prop_str); %grab the time stamp
                                
                                if(ii > 1 && str2num(time) == 0) % if there is no time stamp (i.e == 0)
                                    prop_str = 'KineticCycleTime';%prep cycle time variable
                                    [rc,time_raw] = atsif_getpropertyvalue(signal,prop_str); %grab cycle time variable
                                    time = (str2num(time_raw) * (ii-1)) * 1e6; % calculate time based on cycle time (in microseconds)
                                    frame_time(1) = time; % store calculated time in frame_time array
                                else
                                    frame_time(1) = str2num(time); % store the measured time vector in frame_time array
                                end
                                
                                if(ii == 1) % get the true time for the first frame
                                    [rc,format_time] = atsif_getpropertyvalue(signal,'FormattedTime');
                                    frame_info{1,1} = 'Formatted Time';
                                    frame_info{1,2} = format_time;
                                    time_zero = datenum(format_time,'ddd mmm dd HH:MM:SS yyyy');
                                end
                                if(mod(ii,100) == 0)
                                    fprintf('%i\n',ii)
                                end
                                frame_time(2) = time_zero + frame_time(1) / (24 * 60 * 60 * 1e6) - (10/24);% calculate the matlab time
                            end
                        end
                    end
                end
                atsif_closefile;
            else
                disp('Could not load file.  ERROR - ');
                disp(rc);
            end
        end
    end
end

