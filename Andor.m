classdef Andor < Image
    properties
        timeZero        % start of data set (may not need to keep this)
    end
    
    methods
        function [obj] = Andor(andor_file,memoryStep,startFrame,endFrame,type) % constructor (fills in field values with default settings)
            
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                memoryStep = 0;
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
            
            delta = 30; % sets the frame size for fitting (2*delta)
            obj.memoryStep = memoryStep;
            obj.timeUnits = 1e-6; %recorded in microseconds
            obj.filename = andor_file;% stores the filename with the object
            
            %First get data needed on entire dataSet before reading in each
            %frame individually
            [signal,obj.noFrames,and_size,width,height] = Andor.getDataSetInfo(andor_file); %check set exists, frame dimensions, total number of frames)
            obj.dimensions = [height,width];  % assign dimensions to object (used for fitting and plotting)
            [obj.timeZero] = Andor.getTimeZero (signal); % get the initial time (used to build computer clock time vector)
            
            jj = 0;% counter for tracking stored frames
            
            if endFrame > obj.noFrames || endFrame == 0
                endFrame = obj.noFrames;
            end
            
            for ii = startFrame:endFrame
                
                %read each frame individually
                %=============================================%
                %purpose: grabs a single frame and time stamp
                %inputs: flag, size of linearized frame, frame width, frame height, frame number
                %outputs: frame, time stamps
                [imageData,obj.time(ii,:)] = Andor.getFrameInfo(signal,and_size,width,height,ii);
               
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
            
            obj.time(:,2) = obj.time(:,2)+ obj.timeZero; % add initial time to matlab time to get the computer time
            
        end
        %=============================================%
        % methods that bundle analysis methods %
        function [obj] = analyzeAndorData(obj)
            [obj] = calcMean(obj);
            [obj] = calcRMS(obj);
            [obj] = calcRange(obj);
            [obj] = calcFrameRate(obj);
            [obj] = calcFFT(obj);
            [obj] = calcPSD(obj);
        end
        function [obj] = quickAnalyzeAndorData(obj)
        end        
        %=============================================%
        % methods that plot data from object or find frames which can be plotted%
        function [obj] = inspectFrame(obj,frameNumber)
            [signal,~,and_size,width,height] = Andor.getDataSetInfo(obj.filename);
            [imageData,~] = Andor.getFrameInfo(signal,and_size,width,height,frameNumber);
            obj.tempFrame = imageData;
            figure
            imagesc(obj.tempFrame)
        end
    end
    
    methods(Static)
        %=============================================%
        % functions that recover the image and time stamps %
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
        %=============================================%
        % functions that condition frame for reading in and trimming%
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
                        frame_info = cell(1,2);
                    end
                end
            end
        end
        %=============================================%
        %old function for reading andor data (not used)%
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

