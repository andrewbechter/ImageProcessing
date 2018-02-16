classdef Andor
    properties
        filename        % filename for later reference
        timeZero        % start of data set (may not need to keep this)
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
        
    end
    
    methods
        function [obj] = Andor(andor_file,memoryStep,endFrame,type) % constructor (fills in field values with default settings)
            
            % check the number of input variables and handle missing values
            if nargin <1
                disp('you must specify a data path')
            elseif nargin < 2
                memoryStep = 0;
                type = 'fast';
                endFrame = 0;
            elseif nargin < 3
                type = 'fast';
                endFrame = 0;
            elseif nargin < 4
                type = 'fast';
            end
            
            delta = 50; % sets the frame size for fitting (delta*PSF size)
            obj.memoryStep = memoryStep;
            obj.timeUnits = 1e-6; %recorded in microseconds
            obj.filename = andor_file;% stores the filename with the object
            
            %First get data needed on entire dataSet before reading in each
            %frame individually
            [signal,obj.noFrames,and_size,width,height] = Andor.getDataSetInfo(andor_file); %check set exists, frame dimensions, total number of frames)
            obj.dimensions = [width,height];  % assign dimensions to object (used for fitting and plotting)
            [obj.timeZero] = Andor.getTimeZero (signal); % get the initial time (used to build computer clock time vector)
            
            jj = 0;% counter for tracking stored frames
            
            if endFrame > obj.noFrames || endFrame == 0
                endFrame = obj.noFrames;
            end
            
            for ii = 1:endFrame
                
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
                [obj.iVals(ii,:)] = Andor.xcorrFit(imageData);
                %[iVals(ii,:)] = Andor.coarsefit(imageData);
                
                if strcmp(type,'full')==1
                %=============================================%
                %purpose: check boundaries of cuts vs. the edges of frame
                %inputs: ycenter, xcenter, psf sigma, dimensions, delta
                %outputs: cut locations for frame
                [cuts(ii,:)] = Andor.findFrameCuts(obj.iVals(ii,4),obj.iVals(ii,2),obj.dimensions,delta);
                
                
                %=============================================%
                %inputs: ycenter, xcenter, constant sigma, dimensions, delta
                %outputs: cut locations for frame
                [cutFrame,xdata,flag] = Andor.trimFrame(imageData,cuts(ii,:));%trim the frame
                
                %=============================================%
                %inputs: ycenter, xcenter, constant sigma, dimensions, delta
                %outputs: cut locations for frame
                [obj.fitParameters(ii,:)] = Andor.subPixelPeakV3(cutFrame,obj.iVals(ii,:),0,xdata);% sub pixel fitting
                
               
                
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
        % methods that calculate time series data from object %
        function [obj] = calcFrameRate(obj)
        obj.frameRate = mean(1./diff(obj.timeUnits.*obj.time(:,1)));
        end   
        function [obj] = calcFFT(obj,xdata)
            if nargin <2
               xdata =  obj.fitParameters(:,2);
            end
            Fs = obj.frameRate;
            X = xdata;
            Y = fft(X);
            L = length(X);
            % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            
            % Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
            f = Fs*(0:(L/2))/L;
            obj.FFT(:,1) = f;
            obj.FFT(:,2) = P1;
            
        end 
        function [obj] = calcPSD(obj,xdata)
            if nargin <2
                xdata =  obj.fitParameters(:,2);
            end
            Fs = obj.frameRate;
            N = length(xdata);
            xdft = fft(xdata);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            f = Fs*(0:(N/2))/N;
            obj.PSD(:,1) = f;
            obj.PSD(:,2) = psdx;
        end
        
        %=============================================%
        % methods that calculate image quality data from object %
        function [obj] = calcWFE(obj)
            
        end
        function [obj] = calcStrehlRatio(obj)
        end
        function [obj] = calcFiberCoupling(obj)
        end
        %=============================================%
        % methods that calculate statistics data from object %
        function [obj] = calcMean (obj)
            obj.Mean = mean(obj.fitParameters); % Mean value of Centroid parameters
        end
        function [obj] = calcRange (obj)
            obj.Range = mean(obj.fitParameters); % Mean value of Centroid parameters
        end
        function [obj] = calcRMS(obj)
            obj.RMS = std(obj.fitParameters); % RMS value of Centroid parameters
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
        % methods that plot data from object %
        function psfPlot(obj,data_number)
            InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
            if nargin <2
                data_number = 1;
            end
            ii= obj.storedNums(data_number);% assign frame number of stored frame
            %condition variables for plotting fit slices
            data = obj.frame(:,:,data_number);
            x(1) = obj.fitParameters(ii,1); %peak
            x(2) = obj.fitParameters(ii,2); %xcen
            x(3) = obj.fitParameters(ii,3); %xwidth
            x(4) = obj.fitParameters(ii,4); %ycen
            x(5) = obj.fitParameters(ii,5); %ywidth
            x(6) = obj.fitParameters(ii,6); %rotation
            x(7) = obj.fitParameters(ii,7); %constant
            dispangle = x(6)*180/pi; % convert rotation angle to degrees
            
            if isempty(obj.cuts) == 1
                [X,Y] = meshgrid(1:size(data,2),1:size(data,1));
                Z = data;
            else
                cuts = obj.cuts(data_number,:);
                [X,Y] = meshgrid(cuts(1):cuts(2),cuts(3):cuts(4));
                Z = data(cuts(3):cuts(4),cuts(1):cuts(2)); % cut data at the locations corresponding to the chosen frame size.
            end
            %figure creation
            
            hf2 = figure(ii);
            set(hf2, 'Position', [20 20 950 900])
            alpha(0)
            subplot(4,4,[5,6,7,9,10,11,13,14,15])
            imagesc(X(1,:),Y(:,1)',Z)
            axis equal
            set(gca,'YDir','normal')
            xlabel('X-Pixel')
            ylabel('Y-Pixel')
            set(gca,'Fontsize',16)
            colormap('bone')
            
            % -----Calculate cross sections-------------
            %generate points along horizontal axis
            m = -tan(x(6));% Point slope formula
            b = (-m*x(2) + x(4));
            xvh = X(1,:);
            yvh = xvh*m + b;
            hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
            %generate points along vertical axis
            mrot = -m;
            brot = (mrot*x(4) - x(2));
            yvv = Y(:,1)';
            xvv = yvv*mrot - brot;
            vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);
            
            hold on % Indicate major and minor axis on plot
            %plot points
            plot(xvh,yvh,'r.')
            plot(xvv,yvv,'g.')
            %plot lines
            plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r')
            plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'color',[0 .4 0])
            
            hold off
            
            %fit slices
            ymin = 0;
            ymax = x(1);
            
            xdatafit = linspace(min(X(1,:)), max(X(1,:)),500);
            hdatafit = x(7)+ x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
            
            ydatafit = linspace(min(Y(:,1)), max(Y(:,1)),500);
            vdatafit = x(7)+x(1)*exp(-(ydatafit-x(4)).^2/(2*x(5)^2));
            
            subplot(4,4,[1:3])
            xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
            plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
            axis([min(X(1,:)), max(X(1,:)), ymin, ymax*1.5])
            set(gca,'Fontsize',16)
            
            subplot(4,4,[8,12,16])
            xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
            plot(vPoints,xposv,'b.',vdatafit,ydatafit,'black')
            axis([ymin*1.1 ymax*1.5 min(Y(:,1)), max(Y(:,1))])
            figure(gcf) % bring current figure to front
            set(gca,'Fontsize',16)
            hold off
            
            subplot(4,4,4)
            text(0,+0.5,{['Amplitude:',num2str(x(1),'%100.3f')],['X-center:',num2str(x(2),'%100.3f')],...
                ['\sigma_x (pix):',num2str(x(3),'%100.3f')],['Y-center:',num2str(x(4),'%100.3f')],...
                ['\sigma_y (pix):',num2str(x(5),'%100.3f')],['Angle (deg):',num2str(dispangle,'%100.3f')]...
                ,['Constant:',num2str(x(7),'%100.3f')]},'FontSize',12);
            axis off
        
        end
        function histPlot(obj,x,y,stats)

            if nargin < 2
                x = obj.fitParameters(:,2); %xcen
                xmu = obj.Mean(1,2); %x average
                xrms = obj.RMS(1,2); %x rms
                xrng = obj.Range(1,2); %x range
                y = obj.fitParameters(:,4); %ycen
                ymu = obj.Mean(1,4); %x average
                yrms = obj.RMS(1,4); %x rms
                yrng = obj.Range(1,4); %y range
                n = size(obj.fitParameters(:,1),1); %number of datapoints
            else
                xmu = stats(1); %x average
                xrms = stats(2); %x rms
                xrng = stats(3); %x range
                ymu = stats(4); %x average
                yrms = stats(5); %x rms
                yrng = stats(1); %y range
                n = size(x,1); %number of datapoints
            end
            

            %figure creation
            hf2 = figure;
            set(hf2, 'Position', [20 20 950 950])
            alpha(0)
            p1 = subplot(4,4,[5,6,7,9,10,11,13,14,15]);
            plot(x,y,'.')
            xlim([xmu-4*xrms, xmu+4*xrms])
            ylim([ymu-4*yrms, ymu+4*yrms])
            xlabel('Horizontal Pixels')
            ylabel('Vertical Pixels')
            axis equal
            box on
            p1.FontSize = 12;
            hold on
            
            %Histogram plots
            p2 = subplot(4,4,(1:3));
            histogram(x,round(sqrt(n)))
            xlim([xmu-4*xrms, xmu+4*xrms])
            p2.FontSize = 12;
            
            p3 = subplot(4,4,[8,12,16]);
            histogram(gca,y,round(sqrt(n)))
            xlim([ymu-4*yrms, ymu+4*yrms])
            p3.FontSize = 12;
            set(gca,'view',[90 -90])
            figure(gcf) % bring current figure to front
            hold off
            
            subplot(4,4,4)
            text(0,+0.5,{['Scatter Statistics:'],['\sigma_x: ',num2str(xrms,'%100.2f'),' pix'],['\sigma_y: ',num2str(yrms,'%100.2f'),' pix'],...
                ['Range x: ',num2str(xrng,'%100.2f'),' pix'],['Range y: ',num2str(yrng,'%100.2f'),' pix']},'FontSize',14);
            axis off
            
        end
        function FFTPlot(obj)
            hf = figure();
            set(hf, 'Position', [20 20 1300 350])
            plot(obj.FFT(:,1),obj.FFT(:,2).^2)
            xlabel('f (Hz)')
            ylabel('|P1(f)| (pix)')
            set(gca,'FontSize',16)
            xlim([0.5 obj.frameRate/2])
        end
        function PSDPlot(obj)
            hf = figure();
            set(hf, 'Position', [20 20 1300 350])
            plot(obj.PSD(:,1),obj.PSD(:,2).^2)
            xlabel('f (Hz)')
            ylabel('|P1(f)|^2/f (pix^2/Hz)')
            set(gca,'FontSize',16)
            xlim([0.5 obj.frameRate/2])
        end
        
    end
    
    methods(Static)
        %=============================================%
        % functions that fit image data %
        function [initialValues] = xcorrFit(frame)
            
            s = size(frame);
            % make template image
            center = [s(2) s(1)]/2;
            sig = 5; 
            [MatX,MatY]=meshgrid(1:s(2),1:s(1));
            template=circ_gauss(MatX,MatY,sig,center);
            
%             C = xcorr2(template,frame);
            C = xcorr2_fft(template,frame);
            
            [~,I] = max(C(:)); %find max of linearized data in xcorr
            [cy0,cx0] = ind2sub(size(C),I); % find x,y index locations of maximum count using linear index
            
            shift = [cy0 cx0]-(s-1); % shift( rows (y), columns (x)) +ve is up and left with 0,0 in bottom left
            centroid = round(s/2)-shift; %(y,x)
            
            amp = frame(centroid(1),centroid(2));
            
            rstartpoints = [amp centroid(2) sig 0];
            cstartpoints = [amp centroid(1) sig 0];
            
            geqn = 'a.*exp(-0.5*((x-b)/c).^2)+d';
            
            lb = [0,0,0,0];
            ub = [inf,inf,inf,inf]; % force sigma to be small
            optionsR = fitoptions(geqn);
            optionsR.StartPoint = rstartpoints;
            optionsR.Lower = lb;
            optionsR.Upper = ub;

            optionsC= fitoptions(geqn);
            optionsC.StartPoint = cstartpoints;
            optionsC.Lower = lb;
            optionsC.Upper = ub;
        
            cenrow = frame(centroid(1),:);
            cencol = frame(:,centroid(2));
%             
%             [rfit,rgof] = fit((1:s(2))',cenrow',geqn,optionsR);
%             [cfit,cgof] = fit((1:s(1))',cencol,geqn,optionsR);
            
            [rfit,rgof] = fit((centroid(2)-sig:centroid(2)+sig)',cenrow([centroid(2)-sig:centroid(2)+sig])',geqn,optionsR);
            [cfit,cgof] = fit((centroid(1)-sig:centroid(1)+sig)',cencol([centroid(1)-sig:centroid(1)+sig]),geqn,optionsR);
            
            c=coeffvalues(cfit);
            r=coeffvalues(rfit);
            
            sigma_guess =([r(3),c(3)]);
            initialValues = [amp,centroid(2),sigma_guess(1),centroid(1),sigma_guess(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
            
            
        end
        function [initialValues] = coarsefit(frame)
            
            rsum = sum(frame,1)';
            csum = sum(frame,2);
            
            geqn = 'a.*exp(-0.5*((x-b)/c).^2)+d';
            
            [rmax,rloc] = max(rsum);
            [cmax,cloc] = max(csum);
            
            rstartpoints = [rmax rloc 10 0];
            cstartpoints = [cmax cloc 10 0];
            
            peakamp = max(max(frame)); % find maximum count pixel use amplitude as fit parameter
            [~,I] = max(frame(:)); %find max of linearized data
            [I1,I2] = ind2sub(size(frame),I); % find x,y index locations of maximum count using linear index
            x_coord = I2; %guess for x_coord centroid
            y_coord = I1; %guess for y_coord centroid
            
            rstartpoints = [peakamp x_coord 1 0];
            cstartpoints = [peakamp y_coord 1 0];
            
            lb = [0,0,0.5,0];
            ub = [inf,inf,10,inf]; % force sigma to be small
            optionsR = fitoptions(geqn);
            optionsR.StartPoint = rstartpoints;
            optionsR.Lower = lb;
            optionsR.Upper = ub;
            
            optionsC= fitoptions(geqn);
            optionsC.StartPoint = cstartpoints;
            optionsC.Lower = lb;
            optionsC.Upper = ub;
            
            [rfit,rgof] = fit((1:length(rsum))',rsum,geqn,optionsR);
            [cfit,cgof] = fit((1:length(csum))',csum,geqn,optionsR);
            
            c=coeffvalues(cfit);
            r=coeffvalues(rfit);
            
            center_guess = [cloc,rloc];
            xcen = cloc;
            ycen = rloc;
            sigma_guess =([c(3),r(3)]);
            amp_guess = frame(cloc,rloc);
            initialValues = [amp_guess,xcen,sigma_guess(1),ycen,sigma_guess(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
            
        end
        function [cuts] = findFrameCuts(rcenter,ccenter,dimensions,delta)
            left = max(floor(rcenter-delta),1);
            right = min(floor(rcenter+delta),dimensions(1,2));
            top = max(floor(ccenter-delta),1);
            bot = min(floor(ccenter+delta),dimensions(1,1));
            cuts = [left,right,top,bot];
        end
        function [x] = subPixelPeakV3(frame,x0,FitForOrientation,xdata)
            
            %LSQ fitting for 2D fits
            %=============================================%
            %input format x0 = [Amp,xo,wx,yo,wy,theta,offset];
            options=optimset('Diagnostics','off','Display','none'); % options set for LSQ fit
            
            if FitForOrientation == 1
                lb = [0,0,0,0,0,-pi/4,0];
                ub = [x0(1)*1.2,inf,inf,inf,inf,pi/4,inf];
                % call LSQ curve fit with theta
                [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,frame,lb,ub,options);
            else
                lb = [0,0,0,0,0,0,0];
                ub = [x0(1)*1.2,inf,inf,inf,inf,0,inf];
                % call LSQ curve fit without theta
                [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunction,x0,xdata,frame,lb,ub,options);
                
            end
        end
        %=============================================%
        % functions that recover the time stamps %
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
                imageData = flipud(rot90(imageData));
                
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
        % functions that condition frame %
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
        function [frame,xdata,flag] = trimFrame(frame,cuts)
            [X,Y] = meshgrid(cuts(1):cuts(2),cuts(3):cuts(4));
            % Meshgrid steps over '0' which makes the frame 1 pixel larger than
            % desired. Truncates to correct size. Offsets all frame values (i.e pixels)
            % to be centered at location of PSF - represents the actual detector location
            xdata(:,:,1) = X; %layer 1 is X search space
            xdata(:,:,2) = Y; %layer 2 is Y search space
            %
            Top = Y(1,1);
            Bot = Y(end,1);
            Left = X(1,1);
            Right = X(1,end);
            
            %%%%%%%%%%%----------------------remove bad frames----------------------%%%%%%%%%%%
            flag = 0; % error check
            if(Bot < 1 || Top > size(frame,1) || Left < 1 || Right > size(frame,2))
                flag = 1;
            else
                frame = frame(cuts(3):cuts(4),cuts(1):cuts(2)); % cut data at the locations corresponding to the chosen frame size.
            end
        end
        %=============================================%
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

function F = D2GaussFunction(x,xdata)
% x(1) = 1/sqrt(2*pi*x(3)*x(4));
F = x(7)+ x(1)*exp( -((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ));
end
function F = D2GaussFunctionRot(x,xdata)

xdatarot(:,:,1)= xdata(:,:,1)*cos(x(6)) - xdata(:,:,2)*sin(x(6));
xdatarot(:,:,2)= xdata(:,:,1)*sin(x(6)) + xdata(:,:,2)*cos(x(6));
x0rot = x(2)*cos(x(6)) - x(4)*sin(x(6));
y0rot = x(2)*sin(x(6)) + x(4)*cos(x(6));

F = x(7)+ x(1)*exp(-((xdatarot(:,:,1)-x0rot).^2/(2*x(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*x(5)^2) )    );
end
function F=circ_gauss(X,Y,Sigma,center)
%--------------------------------------------------------------------------
% circ_gauss function                                                General
% Description: Calculate 2D circular Gaussian in a 2-D grid.
% Input  : - Scalar, vector or matrix of X-coordinates in which to calculate
%            the 2-D Gaussian.
%          - same as the x-ccordinates, but for the y-axis.
%          - Sigma of the Gaussian or [SigmaX, SigmaY] in case sigma
%            is different for each axis.
%            By default SigmaY=SigmaX.
%            If empty matrix use default.Si
%          - Center of the Gaussian [X, Y].
%            By default Y=X.
%            Default is [0 0].
%            If empty matrix use default.
%          - Maximum radius of Gaussian behond to set it to zero.
%            Default is Inf.
%            MaxRad is measured from the center of the kernel and not
%            the center of the Gaussian.


% Example: [MatX,MatY]=meshgrid([-10:1:10],[-10:1:10]);
%          F=circ_gauss(MatX,MatY,[1],[0 0]);
%          surface(F);
%--------------------------------------------------------------------------


SigmaX = Sigma;
SigmaY = SigmaX;

X0 = center(1);
Y0 = center(2);

F =exp(-1./(2.).* ((X-X0).^2./SigmaX.^2 +(Y-Y0).^2./SigmaY.^2));

% 1./(2.*pi.*SigmaX.*SigmaY)

% set elements outside MaxRad to zero:
% if (~isinf(cutoff)),
%    MatR = sqrt(X.^2 + Y.^2);
%    I = find(MatR>cutoff);
%    F(I) = 0;
% end
% 
% if (isnan(Norm)),
%    % do not normalize
% else
%    F = Norm.*F./sumnd(F);
% end
end
function c = xcorr2_fft(a,b)
%XCORR2_FFT Two-dimensional cross-correlation evaluated with FFT algorithm.
%   XCORR2_FFT(A,B) computes the cross-correlation of matrices A and B.
%   XCORR2(A) is the autocorrelation function.
%   
%   When matrices A and B are real, XCORR2_FFT is numerically equivalent to
%   XCORR2 but much faster.
%
%   % Example:
%   a = rand(122); b=rand(332);
%   a = a-mean(a(:));
%   b = b-mean(b(:));
%
%   tic,cl = xcorr2(a,b);toc
%   Elapsed time is 0.223502 seconds.
%   tic,cf = xcorr2_fft(a,b);toc
%   Elapsed time is 0.030935 seconds.
% 
%   max(abs(cf(:)-cl(:)))
%   ans = 4.1922e-13
%
%   Author: Alessandro Masullo, 2015
%   Version 1.2
%
%   See also CONV2, XCORR, XCORR2 and FILTER2.

if nargin == 1
	b = a;
end

% Matrix dimensions
adim = size(a);
bdim = size(b);
% Cross-correlation dimension
cdim = adim+bdim-1;

bpad = zeros(cdim);
apad = zeros(cdim);

apad(1:adim(1),1:adim(2)) = a;
bpad(1:bdim(1),1:bdim(2)) = b(end:-1:1,end:-1:1);
ffta = fft2(apad);
fftb = fft2(bpad);
c = real(ifft2(ffta.*fftb));
end
