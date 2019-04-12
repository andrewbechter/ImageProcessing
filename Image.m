classdef Image
    properties
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
        filteredData    % filtered time series data from FouFilter
        coAdd           % all frames summed
        medback         % median background
        sigback         % background std
        SR              % coAdded Strehl ratio
        SR_blur
        airy            % simulated airy pattern
    end
    
    methods
        
        %%-------------------
        % Compute time series 
        %%-------------------
        
        function [obj] = FrameRate(obj)
            % Computes frame rate from time stamps
            % [image object] = FRAMERATE(image object)
            %
            % Parameters
            % ----------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.timeUnites : float
            %   time units ( microseconds etc...)
            %
            % obj.time : float
            %   time stamps
            %
            %
            % Returns
            % -------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.frameRate : float
            %   approximate frame per second (FPS)
            
            obj.frameRate = mean(1./diff(obj.timeUnits.*obj.time(:,1)));
        end
        
        function [obj] = calcFFT(obj,xdata)
            % Computes FFT on time series data in image object
            % [image object] = calcFFT(image object, time series data)
            %
            % Parameters
            % ----------
            % obj : image object
            %   standard image object created using Image class file
            %
            % xdata : 1 x n array
            %
            % Returns
            % -------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.FFT : n x 2 array
            %   frequency (:,1) and power spectrum (:,2)

            if nargin <2
                xdata =  obj.fitParameters(:,2);  % default time series
                % xdata =  obj.instSR;
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
            
            %%---------------------------
            % Assign parameters to object
            %%---------------------------
            obj.FFT(:,1) = f;
            obj.FFT(:,2) = P1;
            
        end
        
        function [obj] = calcPSD(obj,xdata)
            % Computes PSD on time series data in image object
            % [image object] = calcPSD(image object, FFT data)
            %
            % Parameters
            % ----------
            % obj : image object
            %   standard image object created using Image class file
            %
            % xdata : 1 x n array
            %
            % Returns
            % -------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.PSD : n x 2 array
            %   frequency (:,1) and power spectral density (:,2)
            %see matlab documentation here:
            %https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
            
            if nargin <2
                xdata =  obj.fitParameters(:,2); % default
                %xdata =  obj.instSR;
            end
            
            Fs = obj.frameRate;
            N = length(xdata);
            xdft = fft(xdata);
            xdft = xdft(1:N/2+1);
            
            psdx = (1/(Fs*N)) * abs(xdft).^2; %alternative way (same answer) psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            f = Fs*(0:(N/2))/N; % alternative way (same answer) f2 = 0:Fs/N:Fs/2;
            
            %%---------------------------
            % Assign parameters to object
            %%---------------------------
            
            obj.PSD(:,1) = f;
            obj.PSD(:,2) = psdx;
        end
        
        function [obj] = filter(obj,xdata)%(Outdated!)
            % Creates fouerier filter and filters time series data in image object
            % 
            % [image object] = FILTER(image object, time series data)
            %
            % Parameters
            % ----------
            % obj : image object
            %   standard image object created using Image class file
            %
            % xdata : 1 x n array
            %
            % Returns
            % -------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.filteredData : n x 1 array 
            %   filtered Data 
                        
            % Fourier filter function for time-series signal vector y;
            samplingtime = 1./obj.frameRate; % 'samplingtime' is the total duration of sampled signal in sec;
            centerfrequency = 0; % 'centerfrequency' and 'frequencywidth' are the center frequency and width in Hz
            frequencywidth = 0;
            shape = 0; % 'Shape' determines the sharpness of the cut-off. If shape = 1, the filter is Gaussian; as
            % shape increases the filter shape becomes more and more rectangular.
            mode = 0; % Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.
            [obj.filteredData] = FouFilter(xdata,samplingtime,centerfrequency,frequencywidth,shape,mode);
            
        end
        
        %%-------------------
        % Image quality 
        %%-------------------
        
        function [obj] = calcCircPSF(obj,sigma)
            % Finds all the circular PSFs and returns their index in obj.circPSF.
            
            % Computes frame rate from time stamps
            % [image object] = CALCCIRCPSF(image object, wdith parameter)
            %
            % Parameters
            % ----------
            % obj : image object
            %   standard image object created using Image class file
            %
            % sigma : float
            %   gaussian width parameter
            %
            % Returns
            % -------
            % obj : image object
            %   standard image object created using Image class file
            %
            % obj.circPSF : int
            %   index of cicular PSFs  
            %
            % Notes
            % -------
            % So far a sigma value of <5 seems to correlate with a good 'core'
            % so that is the default setting for demonstrator data. 
            
            if nargin < 2
                value = 5;
            else
                value = sigma;
            end
            [ind] = find(obj.fitParameters(:,3) < value & obj.fitParameters(:,5) < value & obj.fitParameters(:,1) > 4500);
            obj.circPSF = ind;
        end
        
        function [obj] = calcFiberCoupling(obj)
        end
        
        function [obj] = addFrames(obj)
            obj.coAdd = sum(obj.frame,3);
        end
        
        function [obj] = calcStrehl(obj,opts)
            
            lambda = opts.lambda; % wavelegnth in meters
            
            focal_length = opts.focal_length; % effective focal length
            
            d = opts.d; % beam diameter
            
            pix = opts.pix; %pixel size
            
            npup = opts.npup; %number of samples across the pupil
            
            alpha = opts.alpha; %blocking fraction
            
            F = focal_length/d; % compute F/number
            
            q = F*lambda/(pix); %samples across PSF
            
            [airy] = Image.diffraction_pattern(q,lambda,npup,alpha); % compute diffraction pattern
            
            [airy,~,~] = Image.centroid_center(airy, 0.1, 5, false); % center image and calculate dx,dy
            
            obj.SR = Image.calculate_SR (obj.coAdd,airy,true);
            
            %obj.imStack = Image.centroid_center(obj.imStack,0.1,5,false);
            
            %obj.SR_blur = Image.calculate_SR (obj.imStack,airy,false);
            
        end
        
        %%-------------------
        % Statistics 
        %%-------------------
        
        function [obj] = calcMean (obj)
            a = obj.fitParameters(:,1)~=0;
            obj.Mean = mean(obj.fitParameters(a,:),1); % Mean value of Centroid parameters
        end
        
        function [obj] = calcRange (obj)
            a = obj.fitParameters(:,1)~=0;
            obj.Range = abs(max(obj.fitParameters(a,:),[],1)-min(obj.fitParameters(a,:),[],1)); % Mean value of Centroid parameters
        end
        
        function [obj] = calcRMS (obj)
            a = obj.fitParameters(:,1)~=0;
            obj.RMS = std(obj.fitParameters(a,:),1); % RMS value of Centroid parameters
        end
        
        function [obj] = calcDelta(obj)
            x = abs(obj.fitParameters(:,2)-obj.Mean(:,2)); % absolute value from the mean in x
            y = abs(obj.fitParameters(:,4)-obj.Mean(:,4)); % absolute value from the mean in y
            r = sqrt(x.^2+y.^2); % pythag to find radius
            obj.delta = r; % assign value to obj
        end
        
        %%-------------------
        % Plotting
        %%-------------------
        
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
            
            %             if isempty(obj.cuts) == 1
            %                 [X,Y] = meshgrid(1:size(data,2),1:size(data,1));
            %                 Z = data;
            %             else
            %                 cuts = obj.cuts(data_number,:);
            %                 [X,Y] = meshgrid(cuts(1):cuts(2),cuts(3):cuts(4));
            %                 Z = data(cuts(3):cuts(4),cuts(1):cuts(2)); % cut data at the locations corresponding to the chosen frame size.
            %             end
            %figure creation
            
            X = obj.gridX(:,:,data_number);
            Y = obj.gridY(:,:,data_number);
            Z = obj.frame(:,:,data_number);
            
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
        
        function masPlot(obj,x,y,stats)
            
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
            
            a=2.33;
            x = a*x;
            xmu = a*xmu;
            xrms = a*xrms; %x rms
            xrng = a*xrng; %x range
            y = a*y; %ycen
            ymu = a*ymu; %x average
            yrms = a*yrms; %x rms
            yrng = a*yrng; %y range
            
            x = x-xmu;
            y = y-ymu;
            
            %figure creation
            hf2 = figure;
            set(hf2, 'Position', [20 20 950 950])
            alpha(0)
            p1 = subplot(4,4,[5,6,7,9,10,11,13,14,15]);
            plot(x,y,'.','color',[0.3,0.75,0.93])
            %             c = viridis(length(x));
            %             for ii = 1:10:length(x)
            %             plot3(x(ii),y(ii),ii,'.','color',c(ii,:))
            %             hold on
            %             end
            
            xlim([-4*xrms, 4*xrms])
            ylim([-4*yrms, 4*yrms])
            xlabel('Centroid scatter (mas)')
            ylabel('Centroid scatter (mas)')
            box on
            p1.FontSize = 16;
            hold on
            
            %Histogram plots
            p2 = subplot(4,4,(1:3));
            histogram(x,round(sqrt(n)),'FaceColor',[0.3,0.75,0.93])
            xlim([-4*xrms, 4*xrms])
            p2.FontSize = 16;
            
            p3 = subplot(4,4,[8,12,16]);
            histogram(gca,y,round(sqrt(n)),'FaceColor',[0.3,0.75,0.93])
            xlim([-4*yrms, 4*yrms])
            p3.FontSize = 16;
            set(gca,'view',[90 -90])
            figure(gcf) % bring current figure to front
            hold off
            
            subplot(4,4,4)
            text(0,+0.5,{['Scatter Statistics:'],['\sigma_x: ',num2str(xrms,'%100.2f'),' mas'],['\sigma_y: ',num2str(yrms,'%100.2f'),' mas'],...
                ['Range x: ',num2str(xrng,'%100.2f'),' mas'],['Range y: ',num2str(yrng,'%100.2f'),' mas']},'FontSize',14);
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
            hf = figure(2);
            hold on
            set(hf, 'Position', [20 20 1300 350])
            semilogy(obj.PSD(:,1),(obj.PSD(:,2)))
            xlabel('f (Hz)')
            ylabel('Power/ Frequency (Pixels^2/Hz)')
            set(gca,'FontSize',16)
            xlim([0.5 obj.frameRate/2])
        end
        
        function sumPlot(obj)
            figure()
            imagesc(obj.coAdd)
        end
        
        %%-------------------
        % Checks 
        %%-------------------
        
        function [obj] = checkFrame(obj, index)
            %this method is useful for finding examples of frames that
            %meet some condition and are stored (i.e can be viewed). For
            %example the default finds circular PSFs (circPSF) stored in
            %the object and adjusts the index value to work with psfPlot
            %method. i.e. obj.psfPlot(exFrames(2))
            
            %inputs : Andor object, index array (e.g. Aus6, Aus6.circPSF)
            %output : populates exFrame property. (example Frames)
            
            if nargin <2
                index = obj.circPSF;
            end
            
            ind = find(mod(index,obj.memoryStep)==0); % find all frames in circPSF index that are also stored
            obj.exFrames = (index(ind)/obj.memoryStep)+1; % adjust index to be used in psfPlot
        end
      
        %%-------------------
        % Saving 
        %%-------------------

        function saveToStruct(obj, filename)
            % this function saves the object made from Image class as a
            % standard matlab structure
            varname = inputname(1);
            props = properties(obj);
            for p = 1:numel(props)
                s.(props{p})=obj.(props{p});
            end
            eval([varname ' = s'])
            save(filename, varname)
        end
        
    end
    
    methods(Static)
        
        %%-------------------
        % Strehl Ratio 
        %%-------------------
        
        function strehl_ratio = calculate_SR (frame,airy,verbose)
            
            %=================================================
            % Find the center of the simulated PSF
            filt_airy = medfilt2(airy);
            
            [num,idx] = max(filt_airy(:));
            
            [airy_x, airy_y] = ind2sub(size(airy), idx);
            
            
            %=================================================
            if verbose
                fprintf( '===============================');
                fprintf( '\n       Center of unaberated PSF: %d , %d \n', airy_x, airy_y );
            end
            
            %=================================================
            % Take a box around the center of the simulated PSF to match
            % measured data
            
            sub_frame = frame; % pass in a centered frame of the correct size
            
            sub_airy = Image.pad(airy,size(frame,1),size(frame,2));
            
            %%---------------------------
            % Gaussian Fit simulated PSF
            %%---------------------------
            
            fitGrid = 3;
            
            sizex = size(sub_airy,2);
            
            sizey = size(sub_airy,1);
            
            [X,Y] = meshgrid(1:sizex,1:sizey);
            
            [fitframe] = Image.pad(sub_airy,fitGrid,fitGrid);% standard size if frame is larger need offsets in x,y
            
            xdata(:,:,1)= Image.pad(X,fitGrid,fitGrid); %trim grid
            
            xdata(:,:,2)= Image.pad(Y,fitGrid,fitGrid); %trim grid
            
            xcen = X(1,ceil(end/2));
            
            ycen = Y(ceil(end/2), 1);
            
            sig = size(fitframe)./fitGrid;
            
            amp = max(max(fitframe));
            
            iVals = [amp,xcen,sig(1),ycen,sig(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
            
            [params1] = Image.subPixelPeakV3(fitframe,iVals,0,xdata);% sub pixel fitting
            
            amp1 = params1(1);
            
            %=================================================
            % Plot simulated PSF
            if verbose
                figure;
                imagesc((sub_airy)./(sum(sum(sub_airy))));
            end
            
            filt_frame = medfilt2(frame);
            
            [num,idx] = max(filt_frame(:));
            
            [obs_x, obs_y] = ind2sub(size(frame),idx);
            
            %%---------------------------
            % Gaussian Fit measured PSF
            %%---------------------------
            
            sizex = size(sub_frame,2);
            
            sizey = size(sub_frame,1);
            
            [X,Y] = meshgrid(1:sizex,1:sizey);
            
            [fitframe] = Image.pad(frame,fitGrid,fitGrid);% standard size if frame is larger need offsets in x,y
            
            xdata(:,:,1)= Image.pad(X,fitGrid,fitGrid); %trim grid
            
            xdata(:,:,2)= Image.pad(Y,fitGrid,fitGrid); %trim grid
            
            xcen = X(1,ceil(end/2));
            
            ycen = Y(ceil(end/2), 1);
            
            sig = size(fitframe)./5;
            
            amp = max(max(fitframe));
            
            iVals = [amp,xcen,sig(1),ycen,sig(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
            
            [params2] = Image.subPixelPeakV3(fitframe,iVals,0,xdata);% sub pixel fitting
            
            amp2 = params2(1);
            
            if verbose
                %=================================================
                fprintf( '===============================');
                fprintf( '\n       Center of aberrated PSF: %d , %d \n', obs_x, obs_y );
                %=============================================
            end
            
            if verbose
                %sub_frame(sub_frame<0) = 0;
                figure;
                imagesc((sub_frame)./(sum(sum(sub_frame))));
                %                 imagesc(log(abs(sub_frame./(sum(sum(sub_frame))))));
                
            end
            
            
            %=============================================
            % Find the maxima of both
            
            max_unobs = max( max( sub_airy ) ) ;
            
            %max_unobs = amp1;
            
            max_obs = max( max( sub_frame ) ) ;
            
            %max_obs = amp2;
            
            %====================
            % Find the sums
            obs_sum = sum( sum(sub_frame) ) ;
            
            unobs_sum = sum( sum(sub_airy) ) ;
            
            %====================
            % Calculate the Strehl ratio S = ( I(x=0) / SUM(I) ) * ( SUM(P) / P(x=0) )
            strehl_ratio = ( max_obs / obs_sum ) * ( unobs_sum / max_unobs ) ;
            
            if verbose
                %=============================================
                fprintf( '==========================================================');
                fprintf( '\n The calculated Strehl ratio is %f  \n', strehl_ratio) ;
                %=================================================
            end
            %             phi = sqrt(abs(log(strehl_ratio)))*lambda/1064e-9;
            %             strehl_ratio = exp(-(phi^2));
            
        end
        
        function [airy] = diffraction_pattern(q,lambda,npup,alpha)
            
            %---------------------
            % upscale computation
            %---------------------
            
            q = 2*q; % gives higher sampling
            
            npup = 2*npup; % keeps grid size same after binning by factor of 2 (factor of 2 is used to generate upscale simluated nyquist PSFs by minimum amount)
            
            %%-----------------------
            % calculate airy pattern
            %%-----------------------
            
            pup_out = circ(npup,npup/q/2);
            
            pup_in = circ(npup,alpha*npup/q/2);
            
            pupil = pup_out-pup_in;
            
            opd = ones(size(pupil));
            
            P = pupil .* exp(1i * 2*pi/lambda * opd);
            
            psf = fftshift(fft2(P));
            
            psf = psf .* conj(psf);
            
            airy = psf ./ sum(psf(:));
            
            
            %%----------
            % Bin down
            %%----------
            
            p = 2; % default by 2
            
            q = 2; % default by 2
            
            [m,n]=size(airy); %M is the original matrix
            
            airy=sum(reshape(airy,p,[]) ,1 );
            
            airy=reshape(airy,m/p,[]).'; %Note transpose
            
            airy=sum(reshape(airy,q,[]) ,1);
            
            airy=reshape(airy,n/q,[]).'; %Note transpose
            
        end
        
        %%-------------------
        % Calibrate
        %%-------------------
        
        function [frame] = subDark(frame,dark)
            frame = frame - dark;
        end
        
        function [frame,medback,bsig] = subBackground(frame,background)
            if nargin<2
                background = false;
            end
            
            if background == false
                
                %%-----------------------------------
                % Remove background (aperture method)
                %%-----------------------------------
                
                X  = size(frame,1);
                Y  = size(frame,2);
                
                [cols , rows] = meshgrid(1:X,1:Y);
                
                centerX = X/2;
                centerY = Y/2;
                
                innerRadius = X/2.4;
                outerRadius = X/2;
                
                array2D = (rows-centerY).^2 +(cols-centerX).^2;
                
                ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;
                
                medback = median(frame(ringPixels));
                
                bsig = std(frame(ringPixels));
                
                frame = frame-medback;
                
            else
                
                %%------------------------------------
                % Remove background (calibration frame)
                %%------------------------------------
                
                frame = frame - background;
                
                medback = median(background);
                
                bsig = std(background);
                
            end
        end
        
        function [frame] = divFlat(frame)
            % no flat field setup in the lab.
        end
        
        function [frame] = sigmaClipping(frame,background)
            sigma = std(background(:));
            ind =  abs(background-mean(mean(background)))> 3*sigma;
            frame(ind) = mean(mean(frame));
        end
        
        %%-------------------
        % Fitting PSFs 
        %%-------------------
        
        function [initialValues,totalCounts,flag] = xcorrFit(frame)
            
            s = size(frame); %find frame dimensions
            totalCounts = sum(sum(frame));
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
            
            if centroid(1)<5 || centroid(1)> s(1)-5 ||centroid(2)<5 || centroid(2)>s(2)-5
                initialValues = [0,0,0,0,0,0,0];
                flag = 1;
                
            else
                initialValues = [amp,centroid(2),0,centroid(1),0,0,0]; %amp,x,y,sigmax,sigmay,theta,offset
                flag = 0;
            end
            
        end
        
        function [flag,sigma]=fit1D(amp,centroid,frame)
            
            sig = 2;
            
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
            
            cenrow = frame(centroid(2),:);
            cencol = frame(:,centroid(1));
            %
            %             [rfit,rgof] = fit((1:s(2))',cenrow',geqn,optionsR);
            %             [cfit,cgof] = fit((1:s(1))',cencol,geqn,optionsR);
            
            [rfit,rgof] = fit((centroid(2)-sig:centroid(2)+sig)',cenrow([centroid(2)-sig:centroid(2)+sig])',geqn,optionsR);
            [cfit,cgof] = fit((centroid(1)-sig:centroid(1)+sig)',cencol([centroid(1)-sig:centroid(1)+sig]),geqn,optionsR);
            
            c=coeffvalues(cfit);
            r=coeffvalues(rfit);
            %central row fit has x values central column has y values
            
            sigma =([r(3),c(3)]); %x,y
            % initialValues = [amp,centroid(2),sigma_guess(1),centroid(1),sigma_guess(2),0,0]; %amp,x,y,sigmax,sigmay,theta,offset
            flag = 0;
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
            top = max(floor(rcenter-delta),1);
            bot = min(floor(rcenter+delta),dimensions(1,1));
            left = max(floor(ccenter-delta),1);
            right = min(floor(ccenter+delta),dimensions(1,2));
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
                ub = [inf,inf,inf,inf,inf,0,inf];
                % call LSQ curve fit without theta
                [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunction,x0,xdata,frame,lb,ub,options);
                
            end
        end
        
        function [centered_frame, dx, dy] = centroid_center(frame, t, n, verbose)
            % Iteratively centeres PSF in image grid using fourier shifting
            % and flux weighted centroid method
            %
            % [centered frame, dx, dy] = CENTROID_CENTER(image frame, threshold, iterations, verbose)
            %
            % Parameters
            % ----------
            % img :PSF image
            %
            % threshold : float
            %   Clipping threshold relative to peak
            %
            % niter : int
            %   Maximum number of iterations to use
            %
            % verbose : bool, optional
            %   If true, print details to the console. Default is false.
            %
            % Returns
            % -------
            % img_centered : m x n array
            %   Centered and unit normalized PSF
            %
            % dx : float
            %   Measured x-offset (column offset)
            %
            % dy : float
            %   Measured y-offset (row offset)
            
            narginchk(3,4)
            if nargin < 4
                verbose = false;
            end
            
            [nrows, ncols] = size(frame);
            
            xc = floor(ncols/2) + 1;
            yc = floor(nrows/2) + 1;
            
            centered_frame = frame;
            
            dxh=[];
            dyh=[];
            
            for kk = 1:n
                [xo, yo] = Image.centroid(centered_frame.*(centered_frame>max(centered_frame(:)*t)));
                
                dx = xo-xc;
                dy = yo-yc;
                
                dxh = [dxh dx];
                dyh = [dyh dy];
                
                centered_frame = image_shift(centered_frame, -dx, -dy);
                
                dr = sqrt(dx.^2 + dy.^2);
                
                vprintf(verbose, ' Iter=%d : dx=%5.2f dy=%5.2f (delta_iter=%6.3f)\n',kk,sum(dxh),sum(dyh),dr)
                
                if dr<0.001
                    break;
                end
            end
            
            dx = sum(dxh);
            dy = sum(dyh);
            vprintf(verbose, ' PSF centered to %6.4f pixels.\n',dr)
            vprintf(verbose, '------------------------------------------------------------\n ')
            
        end
        
        function [xc, yc] = centroid(frame)
            % CENTROID computes image centroid location using flux weight
            % of entire frame
            % [xc, yc] = CENTROID(img)
            
            [rows, cols] = size(frame);
            
            xc = floor(cols/2) + 1;
            yc = floor(rows/2) + 1;
            
            m  = sum(sum(frame));
            
            if m ~= 0
                mx = sum(frame)  * (1:cols)';
                my = sum(frame') * (1:rows)';
                
                xc = mx / m;
                yc = my / m;
            end
        end
                
        %%-------------------
        % Condition Image 
        %%-------------------
        
        function [frame,xdata,flag] = trimFrame(frame,cuts)
            [X,Y] = meshgrid(cuts(1):cuts(2),cuts(3):cuts(4));
            % Meshgrid steps over '0' which makes the frame 1 pixel larger than
            % desired. Truncates to correct size. Offsets all frame values (i.e pixels)
            % to be centered at location of PSF - represents the actual detector location
            xdata(:,:,1) = X; %layer 1 is X search space
            xdata(:,:,2) = Y; %layer 2 is Y search space
            
            Top = Y(1,1);
            Bot = Y(end,1);
            Left = X(1,1);
            Right = X(1,end);
            
            %%%%%%%%%%%----------------------remove bad frames----------------------%%%%%%%%%%%
            flag = 0; % error check
            if(Bot > size(frame,1) || Top < 1  || Left < 1 || Right > size(frame,2))
                flag = 1;
            else
                frame = frame(cuts(3):cuts(4),cuts(1):cuts(2)); % cut data at the locations corresponding to the chosen frame size.
            end
        end
        
        function [out] = pad(in,npix_rows,npix_cols)
            % PAD zero-pads input matrix to size npix by npix or cuts down to size\
            %
            % out = pad_truncate(in,npix_rows,npix_cols)
            %
            % Parameters
            % ----------
            % in : m x n array
            %   Matrix to be padded
            %
            % npix_rows : int
            %   Desired number of rows in output
            %
            % npix_cols : int, optional
            %   Desired number of columns in output. Default is npix_rows (result will
            %   be square: npix_rows x npix_rows)
            %
            % Returns
            % -------
            % out : npix_rows x npix_cols array
            %   Padded input
            %
            
            narginchk(2,3)
            
            if nargin < 3
                npix_cols = npix_rows;
            end
            
            out = zeros(npix_rows,npix_cols);
            
            [nrows, ncols]  = size(in);
            
            ixc = floor(ncols/2 + 1);
            iyc = floor(nrows/2 + 1);
            
            oxc = floor(npix_cols/2 + 1);
            oyc = floor(npix_rows/2 + 1);
            
            dx = npix_cols-ncols;
            dy = npix_rows-nrows;
            
            if dx<=0
                ix1 = ixc - floor(npix_cols/2);
                ix2 = ix1 + npix_cols - 1;
                ox1 = 1;
                ox2 = npix_cols;
            else
                ix1 = 1;
                ix2 = ncols;
                ox1 = oxc - floor(ncols/2);
                ox2 = ox1 + ncols - 1;
            end
            
            if dy<=0
                iy1 = iyc - floor(npix_rows/2);
                iy2 = iy1 + npix_rows - 1;
                oy1 = 1;
                oy2 = npix_rows;
            else
                iy1 = 1;
                iy2 = nrows;
                oy1 = oyc - floor(nrows/2);
                oy2 = oy1 + nrows - 1;
            end
            out(oy1:oy2, ox1:ox2)  = in(iy1:iy2, ix1:ix2);
        end
        
        %%-------------------
        % Interactive filter 
        %%-------------------
        
        function [ry] = ifilter(ix,iy,icenter,iwidth,ishape,imode,ifilt)
            % Eric push test
            % ifilter(x,y) or ifilter(y) or ifilter([x y]) or
            % ry=ifilter(x,y,center,width,shape,plotmode,filtermode)
            % Keyboard-operated interactive Fourier filter function for
            % time-series signal (x,y), with keyboard controls that
            % allow you to adjust the filter parameters continuously
            % while observing the effect on your signal dynamically. Optinal
            % input arguments set the intital values of center frequency,
            % filter width, shape, plotmode (1=linear; 2=semilog frequency;
            % 3=semilog amplitude; 4=log-log) and filtermode ('Band-pass',
            % 'Lowpass', 'Highpass', 'Bandpass', and 'Band-reject (notch)')
            % 'X' key changes x-axis scale on spectrum plot between frequency
            % and time (period); zeroth harmonic (DC level) skipped in spectrum plot.
            % Returns the filtered signal. Press K to list keyboard commands.h
            % T. C. O'Haver (toh@umd.edu), Version 4, May, 2014.
            %
            % Example 1:
            % Periodic waveform with 3 harmonic components
            % x=0:1000; y=sin(x/10)+sin(x/5)+sin(x);ifilter(x,y);
            %
            % Example 2 uses optional input arguments to set intital values:
            % x=0:(1/8000):.3;
            % y=(1+12/100.*sin(2*47*pi.*x)).*sin(880*pi.*x)+(1+12/100.*sin(2*20*pi.*x)).*sin(2000*pi.*x);
            % ry=ifilter(x,y,132,31,18,3,'Band-pass');
            %
            % KEYBOARD CONTROLS when figure window is topmost:
            % Adjust center frequency.......Coarse: < and >
            %                               Fine: left and right cursor arrows
            % Adjust filter width...........Coarse: / and "
            %                               Fine: up and down cursor arrows
            % Filter shape..................A,Z (A more rectangular, Z more Gaussian)
            % Filter mode...................B=bandpass; N or R=notch (band reject)
            %                               H=High-pass; L=Low-pass
            % Select plot mode..............B=bandpass; N or R=notch (band reject);H=High-pass;
            %                               L=Low-pass; C=Comb pass; V=Comb notch
            % Print keyboard commands.......K  Pints this list
            % Print filter parameters.......Q  Prints input arguments: center,width,shape,plotmode,filtermode
            % Print current settings........T  Prints list of current settings
            % Switch SPECTRUM X-axis scale..X switch between frequency and period x scale on POWER SPECTRA
            % Switch OUTPUT Y-axis scale....Y switch between fixed or variable y scale on output plot
            % Play output as sound..........P or Enter
            % Save output as .mat file......S
            %
            global x y CENTER WIDTH SHAPE PLOTMODE FILTERMODE YMODE XMODE
            format short g
            format compact
            switch nargin
                % 'nargin' is the number of arguments
                case 1
                    datasize=size(ix);
                    if isvector(ix)
                        x=1:length(ix); % Use this only to create an x vector if needed
                        y=ix;
                    else
                        if datasize(1)<datasize(2),ix=ix';end
                        x=ix(:,1);
                        y=ix(:,2);
                    end
                    % Adjust x and y vector shape to 1 x n (rather than n x 1)
                    x=reshape(x,1,length(x));
                    y=reshape(y,1,length(y));
                    % If necessary, flip the data vectors so that x increases
                    if x(1)>x(length(x)),
                        disp('x-axis flipped.')
                        x=fliplr(x);
                        y=fliplr(y);
                    end
                    CENTER=1;
                    WIDTH=10;
                    SHAPE=2;
                    FILTERMODE='Band-pass'; % Starts in band-pass mode
                    PLOTMODE=1;  % Starts in plotmode=1 for linear plot
                case 2
                    x=ix;
                    y=iy;
                    % Set initial values of filter parameters.
                    CENTER=1;
                    WIDTH=10;
                    SHAPE=2;
                    FILTERMODE='Band-pass'; % Starts in band-pass mode
                    PLOTMODE=1;  % Starts in plotmode=1 for linear plot
                case 7
                    x=ix;
                    y=iy;
                    CENTER=icenter;
                    WIDTH=iwidth;
                    SHAPE=ishape;
                    FILTERMODE=ifilt;
                    PLOTMODE=imode;
                otherwise
                    disp('Invalid number of arguments')
            end % switch nargin
            
            ry=ones(size(x));
            YMODE=0;XMODE=0;
            % Adjust x and y vector shape to 1 x n (rather than n x 1)
            x=reshape(x,1,length(x));
            y=reshape(y,1,length(y));
            
            % Plot the signal and its power spectrum
            figure(1)
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
            
            % Attaches KeyPress test function to the figure.
            set(gcf,'KeyPressFcn',@ReadKey)
            uicontrol('Style','text')
            % end of outer function
        end
    end
end

% ----------------------------SUBFUNCTIONS--------------------------------

%%----------------
% FOURIER FILTER 
%%----------------
function ReadKey(obj,eventdata)
% Interprets key presses from the Figure window.
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
global x y CENTER WIDTH SHAPE PLOTMODE FILTERMODE ry YMODE XMODE
key=get(gcf,'CurrentCharacter');
if ischar(key),
    switch double(key),
        case 28
            % Pans "CENTER" one point down when left arrow pressed.
            CENTER=CENTER-1;
            if CENTER<1,CENTER=1;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 29
            % Pans "CENTER" one point up when left arrow pressed.
            CENTER=CENTER+1;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 44
            % Pans "CENTER" 10 points down when < key pressed.
            CENTER=CENTER-10;
            if CENTER<1,CENTER=1;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 46
            % Pans "CENTER" 10 points up when > key pressed.
            CENTER=CENTER+10;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 30
            % Zooms filter width one point up when up arrow pressed.
            WIDTH=WIDTH+1;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 31
            % Zooms filter width one point down when down arrow pressed.
            WIDTH=WIDTH-1;
            if WIDTH<1,WIDTH=1;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 39
            % Zooms filter width up when / pressed.
            WIDTH=WIDTH+10;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 47
            % Zooms filter width down when ' pressed.
            WIDTH=WIDTH-10;
            if WIDTH<1,WIDTH=1;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 122
            % When 'z' key is pressed, shape is made more rectangular
            SHAPE=SHAPE-1;
            if SHAPE<0,SHAPE=0;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 97
            % When 'a' key is pressed, shape is made more Gaussian
            SHAPE=SHAPE+1;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 49 % When '1' key is pressed, set plot mode to linear/linear
            PLOTMODE=1;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 50 % When '2' key is pressed, set plot mode to log x, linear y
            PLOTMODE=2;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 51
            % When '3' key is pressed, set plot mode to linear x, log y
            PLOTMODE=3;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 52
            % When '4' key is pressed,  set plot mode to log y, log x
            PLOTMODE=4;
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 98
            % When 'b' key is pressed, set to bandpass mode
            FILTERMODE='Band-pass';
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case {110, 114}
            % When 'n' or 'r' key is pressed, set to band reject (FILTERMODE) node
            FILTERMODE='Band-reject (notch)';
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 104
            % When 'H' key is pressed, makes high-pass filter
            FILTERMODE='High-pass';
            ry=RedrawFourierFilter(x,y,length(y)/2,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 108
            % When 'L' is pressed, makes low=pass filter
            FILTERMODE='Low-pass';
            ry=RedrawFourierFilter(x,y,0,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 99
            % When 'C' is pressed, makes comb filter
            FILTERMODE='Comb pass';
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 118
            % When 'V' is pressed, makes comb filter
            FILTERMODE='Comb notch';
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 13
            % When 'Enter' is pressed, plays filtered signal as sound
            % without scaliang
            sound(ry,8000);
        case 112
            % When 'p' key is pressed, plays filtered signal as sound, scaled
            % so that the sound is played as loud as possible without clipping.
            soundsc(ry,8000);
        case 115
            % When 's' key is pressed, filtered signal is saved as
            % FilteredOutput.mat
            save FilteredOutput ry
        case 107
            % When 'k' key is pressed, prints out table of keyboard commands
            disp('')
            disp('KEYBOARD CONTROLS when figure window is topmost:')
            disp('Adjust center frequency.......Coarse: < and >')
            disp('                              Fine: left and right cursor arrows')
            disp('Adjust filter width...........Coarse: / and "  ')
            disp('                              Fine: up and down cursor arrows')
            disp('Filter shape..................A,Z (A more rectangular, Z more Gaussian)')
            disp('Filter mode...................B=bandpass; N or R=notch (band reject); H=High-pass;')
            disp('                              L=Low-pass; C=Comb pass; V=Comb notch (reject)')
            disp('Select plot mode..............1=linear; 2=semilog frequency' )
            disp('                              3=semilog amplitude; 4=log-log')
            disp('Print keyboard commands.......K  Pints this list')
            disp('Print filter parameters.......Q  Prints input arguments: center,width,shape,plotmode,filtermode')
            disp('Print current settings........T  Prints list of current settings')
            disp('Switch SPECTRUM X-axis scale..X switch between frequency and period x scale on POWER SPECTRA')
            disp('Switch OUTPUT Y-axis scale....Y switch between fixed or variable y scale on output plot')
            disp('Play output as sound..........P or Enter')
            disp('Save output as .mat file......S')
            disp(' ')
        case 113
            % When 'Q' is pressed, prints fourier filter parameters on a single line
            disp([ num2str(CENTER) ',' num2str(WIDTH)  ',' num2str(SHAPE)  ',' num2str(PLOTMODE) ',' num2str(FILTERMODE) ])
        case 116
            % When 'T' is pressed, prints list of current settings
            disp('------ iFilter Settings --------' )
            disp(['Number of points in signal: ' num2str(length(y)) ] )
            disp(['Duration of signal (x range): ' num2str(range(x)) ] )
            disp(['Interval between x-axis values: ' num2str(x(2)-x(1)) ] )
            disp(['Center harmonic: ' num2str(CENTER) ] )
            disp(['Center frequency: ' num2str(CENTER./range(x)) ] )
            disp(['Center period: ' num2str(range(x)./CENTER) ] )
            disp(['Width in harmonics: ' num2str(WIDTH) ] )
            disp(['Width in frequency: ' num2str(WIDTH./range(x)) ] )
            disp(['Filter mode: ' FILTERMODE ] )
        case 102
            
        case 121
            % When 'Y' is pressed, toggles between fixed or variable y scale on OUTPUT plot
            if YMODE==0,YMODE=1;else YMODE=0;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        case 120
            % When 'X' is pressed, toggles between frequency and period x=vxis on spectrum
            if XMODE==0,XMODE=1;else XMODE=0;end
            ry=RedrawFourierFilter(x,y,CENTER,WIDTH,SHAPE,PLOTMODE,FILTERMODE,YMODE,XMODE);
        otherwise  % if key pressed is unnasigned
            UnassignedKey=double(key) % print unassigned key code
            disp('Press K to print out list of keyboard commands')
    end % switch
end % if
end
function ry=RedrawFourierFilter(xvector,yvector,centerfrequency,filterwidth,filtershape,mode,FILTERMODE,YMODE,XMODE)
% Separate graph windows for the original and filtered signals.
% Computes and plots fourier filter for signal yvector.
% Centerfrequency and filterwidth are the center frequency and
% width of the pass band, in harmonics, 'filtershape' determines
% the sharpness of the cut-off. Plot modes: mode=1 linear x and y;
% mode=2 log x linear y; mode=3 linear x; log y; mode=4 log y log x
fy=fft(yvector);
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape.
if strcmp(FILTERMODE,'Band-pass'),
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'High-pass')
    centerfrequency=length(xvector)/2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Low-pass')
    centerfrequency=0;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Band-reject (notch)')
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    ffilter=1-[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Comb pass')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<30,
        ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
        ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
        n=n+1;
    end
    ffilter=[ffilter1,ffilter2];
end
if strcmp(FILTERMODE,'Comb notch')
    n=2;
    ffilter1=shape(lft1,centerfrequency+1,filterwidth,filtershape);
    ffilter2=shape(lft2,length(fy)-centerfrequency+1,filterwidth,filtershape);
    while n<30,
        ffilter1=ffilter1+shape(lft1,n*(centerfrequency+1),filterwidth,filtershape);
        ffilter2=ffilter2+shape(lft2,length(fy)-n*(centerfrequency+1),filterwidth,filtershape);
        n=n+1;
    end
    ffilter=1-[ffilter1,ffilter2];
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier Transform of signal
ry=real(ifft(ffy));

subplot(3,1,1)  % Plot original signal in top plot
plot(xvector,yvector,'b');
title('iFilter 4:       Press K to display keyboard commands.')
xlabel('INPUT: Original Signal   x=time')
axis([xvector(1) xvector(length(xvector)) min(yvector) max(yvector)]);

subplot(3,1,3)  % Plot filtered signal in lower plot
plot(xvector,ry,'r');
if YMODE,
    axis([xvector(1) xvector(length(xvector)) min(ry) max(ry)]);
else
    axis([xvector(1) xvector(length(xvector)) min(yvector) max(yvector)]);
end
title([ FILTERMODE ' mode:  Freq= ' num2str(centerfrequency./range(xvector)) '    Period= ' num2str(range(xvector)./centerfrequency) '   Width= ' num2str(filterwidth./range(xvector)) '     Shape=  ' num2str(filtershape)])
xlabel('OUTPUT: Filtered signal    x=time')
if YMODE==0,
    ylabel('Fixed Y scale')
else
    ylabel('Variable Y scale')
end
subplot(3,1,2)    % Plot power spectrum and filter in middle plot
py=fy .* conj(fy); % Compute power spectrum
plotrange=2:length(fy)/2;
if XMODE,
    f=range(xvector)./(plotrange-1);
else
    f=((plotrange-1)./range(xvector));
end
switch mode,
    case 1
        plot(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
        ylabel('Linear y')
    case 2
        semilogx(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
        ylabel('Linear y')
    case 3
        semilogy(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
        ylabel('Log y')
    case 4
        loglog(f,real(py(plotrange)),f,max(real(py(plotrange))).*ffilter(plotrange),'r')
        ylabel('Log y')
    otherwise,
end
title('POWER SPECTRA:  BLUE = Input signal    RED = Filter')
if XMODE,
    xlabel('x=Time (e.g. seconds)')
else
    xlabel('x=Frequency (e.g. cycles/second)')
end
% Expand frequency axis if filter is in the lower half of frequency axis
if centerfrequency+filterwidth/2<length(fy)/4,
    axis([0 max(f)/2 min(py) 1.1*max(real(py(plotrange)))])
else
    axis([0 max(f) min(py) 1.1*max(real(py(plotrange)))])
end
end
function g = shape(x,pos,wid,n)
%  shape(x,pos,wid,n) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Lorentzian (1/x^2) when n=0, Gaussian (exp(-x^2))
%  when n=1, and becomes more rectangular as n increases.
%  Example: shape([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
if n==0
    g=ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
else
    g = exp(-((x-pos)./(0.6.*wid)) .^(2*round(n)));
end
end
function r = range(arr)
r = max(arr) - min(arr);
end
function [ry]=FouFilter(y,samplingtime,centerfrequency,frequencywidth,shape,mode)
% Fourier filter function for time-series signal vector y; 'samplingtime'
% is the total duration of sampled signal in sec, millisec, or microsec;
% 'centerfrequency' and 'frequencywidth' are the center frequency and width
% of the filter in Hz, KHz, or MHz, respectively; 'Shape' determines the
% sharpness of the cut-off. If shape = 1, the filter is Gaussian; as
% shape increases the filter shape becomes more and more rectangular.
% Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.
% FouFilter returns the filtered signal.
%
% Example: Sine wave in noisy background.
% First half is just noise; sine wave starts halfway through.
% xx=[0:.001:2*pi]';
% signal=sin(20*xx);
% noise=randn(size(xx));
% x=1:2*length(xx)';
% y=[noise;signal+noise]; % sine wave is added halfway through.
% SignalToNoiseRatio=std(signal)/std(noise)
% FilteredSignal=foufilter(y',1,20,100,5,0);
% subplot(2,1,1)
% plot(x,y);
% title('First half is just noise; sine wave starts halfway through')
% subplot(2,1,2)
% plot(x,FilteredSignal);
% title('Signal filtered with FouFilter.m')
%
%  T. C. O'Haver (toh@umd.edu),  version 1.5, May, 2007

center=centerfrequency*samplingtime; %  center harmonic (fourier component)
width=frequencywidth*samplingtime; %  width of filter (in harmonics)

fy=fft(y); % Fourier transform of signal
lft1=[1:(length(fy)/2)];
lft2=[(length(fy)/2+1):length(fy)];
% Compute filter shape
ffilter1=ngaussian(lft1,center+1,width,shape);
ffilter2=ngaussian(lft2,length(fy)-center+1,width,shape);
ffilter=[ffilter1,ffilter2];
if mode==1,
    ffilter=1-ffilter;
end
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter;  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Recover filter signal from Fourier transform
end

%%---------------
% MATH FUNCTIONS 
%%---------------
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
function F = circ_gauss(X,Y,Sigma,center)
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
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Gaussian when n=1, becomes more rectangular as n increases.
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6005615.*wid)) .^(2*round(n)));
end
function C = circ(N,R)
C = zeros(N);

[X,Y] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);

X = X+0.5; Y=Y+0.5;

Z = sqrt(X.^2 + Y.^2);

C(Z<R) = 1;

end
function [b,bsig] = back(frame)

imageSizeX = size(frame,2); % size of frame cols

imageSizeY = size(frame,1); % size of frame rows

[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the apreture circle in the image.

centerX = imageSizeX/2;

centerY = imageSizeY/2;

innerRadius = imageSizeX/3;

outerRadius = imageSizeX/2;

array2D = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2;

ringPixels = array2D >= innerRadius.^2 & array2D <= outerRadius.^2;
% ringPixels is a 2D "logical" array.

b = mean(frame(ringPixels));

bsig = std(frame(ringPixels));

end

%%-----------------
% FOURIER TRANSFORM
%%-----------------
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
function [psf] = opd2psf(pupil,opd,lambda)

P = pupil .* exp(i * 2*pi/lambda * opd);

psf = fftshift(fft2(P));

psf = psf .* conj(psf);

psf = psf ./ sum(psf(:));

end
function result = image_shift(image,dx,dy,complex)
% IMAGE_SHIFT shifts an image via FFT
%
% result = IMAGE_SHIFT(image,dx,dy,complex)
%
% Shift an image by the specified dx and dy offsets. These offsets may be
% non-integer as the shift is implemented by introducing a Fourier domain
% tilt term.
%
% Parameters
% ----------
% image : m x n array
%   Image to be shifted
%
% dx, dy : float or 1 x n arrays
%   (x, y) or (column, row) shift. If dx and dy are vectors, the result will
%   be a shift-and-add over the set of specified shifts. This is useful for
%   simulating jitter effects or motion blur.
%
% complex : bool, optional
%   If true, the complex result is returned. If false (default), only the
%   real portion of thre result is returned.
%
% Returns
% -------
% result : m x n array
%   Shifted image
%


if nargin < 4
    complex = false;
end

if (dx==0) && (dy==0)
    result = image;
    return;
end

if isempty(dx) || isempty(dy)
    result = image;
    return;
end

[m, n] = size(image);

N = length(dx);

xc = floor(n/2)+1;
yc = floor(m/2)+1;

[X, Y] = meshgrid(1:n,1:m);

X = X - xc;
Y = Y - yc;

%I = fftshift(fft2(image));
I = fft2(image);

T = 0*I;

for k=1:N
    
    px = -2*pi*(X)/n * dx(k);
    py = -2*pi*(Y)/m * dy(k);
    
    T = T + exp(1i * (px + py));
end

%I = I .* T;
I = I .* fftshift(T);

%result = real(ifft2(fftshift(I)));
if complex
    result = ifft2(I);
else
    result = real(ifft2(I));
end

end
function vprintf(verbose, varargin)
% VPRINTF passes fprintf arguments if verbose is true.
%
% vprintf(verbose, <fprintf args>)
%
% Parameters
% ----------
% verbose : bool
%   If true, passes calls fprintf with arguments in varargin
%
% See Also
% --------
% FPRINTF
%


if verbose
    fprintf(varargin{:})
end
end

%%-----------------
% COLOR MAPS
%%-----------------

function C = viridis(N)
%VIRIDIS Blue-green-yellow colour map
%   VIRIDIS(N) returns an N-by-3 matrix containing a colormap.
%   The colors begin with dark purplish-blue and blue, range
%   through green and end with yellow.
%
%   VIRIDIS is the new default colormap for matplotlib
%
% Created by Ed Hawkins (@ed_hawkins) and Kevin Anchukaitis

viridi = [
    0.26700401  0.00487433  0.32941519
    0.26851048  0.00960483  0.33542652
    0.26994384  0.01462494  0.34137895
    0.27130489  0.01994186  0.34726862
    0.27259384  0.02556309  0.35309303
    0.27380934  0.03149748  0.35885256
    0.27495242  0.03775181  0.36454323
    0.27602238  0.04416723  0.37016418
    0.2770184   0.05034437  0.37571452
    0.27794143  0.05632444  0.38119074
    0.27879067  0.06214536  0.38659204
    0.2795655   0.06783587  0.39191723
    0.28026658  0.07341724  0.39716349
    0.28089358  0.07890703  0.40232944
    0.28144581  0.0843197   0.40741404
    0.28192358  0.08966622  0.41241521
    0.28232739  0.09495545  0.41733086
    0.28265633  0.10019576  0.42216032
    0.28291049  0.10539345  0.42690202
    0.28309095  0.11055307  0.43155375
    0.28319704  0.11567966  0.43611482
    0.28322882  0.12077701  0.44058404
    0.28318684  0.12584799  0.44496
    0.283072    0.13089477  0.44924127
    0.28288389  0.13592005  0.45342734
    0.28262297  0.14092556  0.45751726
    0.28229037  0.14591233  0.46150995
    0.28188676  0.15088147  0.46540474
    0.28141228  0.15583425  0.46920128
    0.28086773  0.16077132  0.47289909
    0.28025468  0.16569272  0.47649762
    0.27957399  0.17059884  0.47999675
    0.27882618  0.1754902   0.48339654
    0.27801236  0.18036684  0.48669702
    0.27713437  0.18522836  0.48989831
    0.27619376  0.19007447  0.49300074
    0.27519116  0.1949054   0.49600488
    0.27412802  0.19972086  0.49891131
    0.27300596  0.20452049  0.50172076
    0.27182812  0.20930306  0.50443413
    0.27059473  0.21406899  0.50705243
    0.26930756  0.21881782  0.50957678
    0.26796846  0.22354911  0.5120084
    0.26657984  0.2282621   0.5143487
    0.2651445   0.23295593  0.5165993
    0.2636632   0.23763078  0.51876163
    0.26213801  0.24228619  0.52083736
    0.26057103  0.2469217   0.52282822
    0.25896451  0.25153685  0.52473609
    0.25732244  0.2561304   0.52656332
    0.25564519  0.26070284  0.52831152
    0.25393498  0.26525384  0.52998273
    0.25219404  0.26978306  0.53157905
    0.25042462  0.27429024  0.53310261
    0.24862899  0.27877509  0.53455561
    0.2468114   0.28323662  0.53594093
    0.24497208  0.28767547  0.53726018
    0.24311324  0.29209154  0.53851561
    0.24123708  0.29648471  0.53970946
    0.23934575  0.30085494  0.54084398
    0.23744138  0.30520222  0.5419214
    0.23552606  0.30952657  0.54294396
    0.23360277  0.31382773  0.54391424
    0.2316735   0.3181058   0.54483444
    0.22973926  0.32236127  0.54570633
    0.22780192  0.32659432  0.546532
    0.2258633   0.33080515  0.54731353
    0.22392515  0.334994    0.54805291
    0.22198915  0.33916114  0.54875211
    0.22005691  0.34330688  0.54941304
    0.21812995  0.34743154  0.55003755
    0.21620971  0.35153548  0.55062743
    0.21429757  0.35561907  0.5511844
    0.21239477  0.35968273  0.55171011
    0.2105031   0.36372671  0.55220646
    0.20862342  0.36775151  0.55267486
    0.20675628  0.37175775  0.55311653
    0.20490257  0.37574589  0.55353282
    0.20306309  0.37971644  0.55392505
    0.20123854  0.38366989  0.55429441
    0.1994295   0.38760678  0.55464205
    0.1976365   0.39152762  0.55496905
    0.19585993  0.39543297  0.55527637
    0.19410009  0.39932336  0.55556494
    0.19235719  0.40319934  0.55583559
    0.19063135  0.40706148  0.55608907
    0.18892259  0.41091033  0.55632606
    0.18723083  0.41474645  0.55654717
    0.18555593  0.4185704   0.55675292
    0.18389763  0.42238275  0.55694377
    0.18225561  0.42618405  0.5571201
    0.18062949  0.42997486  0.55728221
    0.17901879  0.43375572  0.55743035
    0.17742298  0.4375272   0.55756466
    0.17584148  0.44128981  0.55768526
    0.17427363  0.4450441   0.55779216
    0.17271876  0.4487906   0.55788532
    0.17117615  0.4525298   0.55796464
    0.16964573  0.45626209  0.55803034
    0.16812641  0.45998802  0.55808199
    0.1666171   0.46370813  0.55811913
    0.16511703  0.4674229   0.55814141
    0.16362543  0.47113278  0.55814842
    0.16214155  0.47483821  0.55813967
    0.16066467  0.47853961  0.55811466
    0.15919413  0.4822374   0.5580728
    0.15772933  0.48593197  0.55801347
    0.15626973  0.4896237   0.557936
    0.15481488  0.49331293  0.55783967
    0.15336445  0.49700003  0.55772371
    0.1519182   0.50068529  0.55758733
    0.15047605  0.50436904  0.55742968
    0.14903918  0.50805136  0.5572505
    0.14760731  0.51173263  0.55704861
    0.14618026  0.51541316  0.55682271
    0.14475863  0.51909319  0.55657181
    0.14334327  0.52277292  0.55629491
    0.14193527  0.52645254  0.55599097
    0.14053599  0.53013219  0.55565893
    0.13914708  0.53381201  0.55529773
    0.13777048  0.53749213  0.55490625
    0.1364085   0.54117264  0.55448339
    0.13506561  0.54485335  0.55402906
    0.13374299  0.54853458  0.55354108
    0.13244401  0.55221637  0.55301828
    0.13117249  0.55589872  0.55245948
    0.1299327   0.55958162  0.55186354
    0.12872938  0.56326503  0.55122927
    0.12756771  0.56694891  0.55055551
    0.12645338  0.57063316  0.5498411
    0.12539383  0.57431754  0.54908564
    0.12439474  0.57800205  0.5482874
    0.12346281  0.58168661  0.54744498
    0.12260562  0.58537105  0.54655722
    0.12183122  0.58905521  0.54562298
    0.12114807  0.59273889  0.54464114
    0.12056501  0.59642187  0.54361058
    0.12009154  0.60010387  0.54253043
    0.11973756  0.60378459  0.54139999
    0.11951163  0.60746388  0.54021751
    0.11942341  0.61114146  0.53898192
    0.11948255  0.61481702  0.53769219
    0.11969858  0.61849025  0.53634733
    0.12008079  0.62216081  0.53494633
    0.12063824  0.62582833  0.53348834
    0.12137972  0.62949242  0.53197275
    0.12231244  0.63315277  0.53039808
    0.12344358  0.63680899  0.52876343
    0.12477953  0.64046069  0.52706792
    0.12632581  0.64410744  0.52531069
    0.12808703  0.64774881  0.52349092
    0.13006688  0.65138436  0.52160791
    0.13226797  0.65501363  0.51966086
    0.13469183  0.65863619  0.5176488
    0.13733921  0.66225157  0.51557101
    0.14020991  0.66585927  0.5134268
    0.14330291  0.66945881  0.51121549
    0.1466164   0.67304968  0.50893644
    0.15014782  0.67663139  0.5065889
    0.15389405  0.68020343  0.50417217
    0.15785146  0.68376525  0.50168574
    0.16201598  0.68731632  0.49912906
    0.1663832   0.69085611  0.49650163
    0.1709484   0.69438405  0.49380294
    0.17570671  0.6978996   0.49103252
    0.18065314  0.70140222  0.48818938
    0.18578266  0.70489133  0.48527326
    0.19109018  0.70836635  0.48228395
    0.19657063  0.71182668  0.47922108
    0.20221902  0.71527175  0.47608431
    0.20803045  0.71870095  0.4728733
    0.21400015  0.72211371  0.46958774
    0.22012381  0.72550945  0.46622638
    0.2263969   0.72888753  0.46278934
    0.23281498  0.73224735  0.45927675
    0.2393739   0.73558828  0.45568838
    0.24606968  0.73890972  0.45202405
    0.25289851  0.74221104  0.44828355
    0.25985676  0.74549162  0.44446673
    0.26694127  0.74875084  0.44057284
    0.27414922  0.75198807  0.4366009
    0.28147681  0.75520266  0.43255207
    0.28892102  0.75839399  0.42842626
    0.29647899  0.76156142  0.42422341
    0.30414796  0.76470433  0.41994346
    0.31192534  0.76782207  0.41558638
    0.3198086   0.77091403  0.41115215
    0.3277958   0.77397953  0.40664011
    0.33588539  0.7770179   0.40204917
    0.34407411  0.78002855  0.39738103
    0.35235985  0.78301086  0.39263579
    0.36074053  0.78596419  0.38781353
    0.3692142   0.78888793  0.38291438
    0.37777892  0.79178146  0.3779385
    0.38643282  0.79464415  0.37288606
    0.39517408  0.79747541  0.36775726
    0.40400101  0.80027461  0.36255223
    0.4129135   0.80304099  0.35726893
    0.42190813  0.80577412  0.35191009
    0.43098317  0.80847343  0.34647607
    0.44013691  0.81113836  0.3409673
    0.44936763  0.81376835  0.33538426
    0.45867362  0.81636288  0.32972749
    0.46805314  0.81892143  0.32399761
    0.47750446  0.82144351  0.31819529
    0.4870258   0.82392862  0.31232133
    0.49661536  0.82637633  0.30637661
    0.5062713   0.82878621  0.30036211
    0.51599182  0.83115784  0.29427888
    0.52577622  0.83349064  0.2881265
    0.5356211   0.83578452  0.28190832
    0.5455244   0.83803918  0.27562602
    0.55548397  0.84025437  0.26928147
    0.5654976   0.8424299   0.26287683
    0.57556297  0.84456561  0.25641457
    0.58567772  0.84666139  0.24989748
    0.59583934  0.84871722  0.24332878
    0.60604528  0.8507331   0.23671214
    0.61629283  0.85270912  0.23005179
    0.62657923  0.85464543  0.22335258
    0.63690157  0.85654226  0.21662012
    0.64725685  0.85839991  0.20986086
    0.65764197  0.86021878  0.20308229
    0.66805369  0.86199932  0.19629307
    0.67848868  0.86374211  0.18950326
    0.68894351  0.86544779  0.18272455
    0.69941463  0.86711711  0.17597055
    0.70989842  0.86875092  0.16925712
    0.72039115  0.87035015  0.16260273
    0.73088902  0.87191584  0.15602894
    0.74138803  0.87344918  0.14956101
    0.75188414  0.87495143  0.14322828
    0.76237342  0.87642392  0.13706449
    0.77285183  0.87786808  0.13110864
    0.78331535  0.87928545  0.12540538
    0.79375994  0.88067763  0.12000532
    0.80418159  0.88204632  0.11496505
    0.81457634  0.88339329  0.11034678
    0.82494028  0.88472036  0.10621724
    0.83526959  0.88602943  0.1026459
    0.84556056  0.88732243  0.09970219
    0.8558096   0.88860134  0.09745186
    0.86601325  0.88986815  0.09595277
    0.87616824  0.89112487  0.09525046
    0.88627146  0.89237353  0.09537439
    0.89632002  0.89361614  0.09633538
    0.90631121  0.89485467  0.09812496
    0.91624212  0.89609127  0.1007168
    0.92610579  0.89732977  0.10407067
    0.93590444  0.8985704   0.10813094
    0.94563626  0.899815    0.11283773
    0.95529972  0.90106534  0.11812832
    0.96489353  0.90232311  0.12394051
    0.97441665  0.90358991  0.13021494
    0.98386829  0.90486726  0.13689671
    0.99324789  0.90615657  0.1439362];

P = size(viridi,1);

if nargin < 1
    N = P;
end

%N = min(N,P);
C = interp1(1:P, viridi, linspace(1,P,N), 'linear');
end