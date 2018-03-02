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
    end
    
    methods
        %=============================================%
        % methods that calculate things using time series data from object %
        function [obj] = calcFrameRate(obj)
            obj.frameRate = mean(1./diff(obj.timeUnits.*obj.time(:,1)));
        end
        function [obj] = calcFFT(obj,xdata)
            if nargin <2
                xdata =  obj.fitParameters(:,4);
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
            %see matlab documentation here:
            %https://www.mathworks.com/help/signal/ug/power-spectral-density-estimates-using-fft.html
            
            if nargin <2
                xdata =  obj.fitParameters(:,4);
            end
            Fs = obj.frameRate;
            N = length(xdata);
            xdft = fft(xdata);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * xdft.*conj(xdft); %alternative way (same answer) psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            f = Fs*(0:(N/2))/N; % alternative way (same answer) f2 = 0:Fs/N:Fs/2;
            obj.PSD(:,1) = f;
            obj.PSD(:,2) = psdx;
        end
        function [obj] = filter(obj,xdata)
            % Fourier filter function for time-series signal vector y;
            samplingtime = 1./obj.frameRate; % 'samplingtime' is the total duration of sampled signal in sec;
            centerfrequency = 0; % 'centerfrequency' and 'frequencywidth' are the center frequency and width in Hz
            frequencywidth = 0;
            shape = 0; % 'Shape' determines the sharpness of the cut-off. If shape = 1, the filter is Gaussian; as
            % shape increases the filter shape becomes more and more rectangular.
            mode = 0; % Set mode = 0 for band-pass filter, mode = 1 for band-reject (notch) filter.
            [obj.filteredData] = FouFilter(xdata,samplingtime,centerfrequency,frequencywidth,shape,mode);
            
        end
        
        
        %=============================================%
        % methods that calculate image quality data from object %
        function [obj] = calcCircPSF(obj,sigma)
            % this function finds all the circular PSFs and returns their index in obj.circPSF.
            % So far a sigma value of <5 seems to correlate with a good 'core'
            % so that is the default setting.
            if nargin < 2
                value = 5;
            else
                value = sigma;
            end
            [ind] = find(obj.fitParameters(:,3) < value & obj.fitParameters(:,5) < value & obj.fitParameters(:,1) > 4500);
            obj.circPSF = ind;
        end
        function [obj] = calcWFE(obj)
        end
        function [obj] = calcStrehlRatio(obj)
        end
        function [obj] = calcFiberCoupling(obj)
        end
        %=============================================%
        % methods that calculate statistics data from object %
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
            obj.RMS = std(obj.fitParameters(a,:),0,1); % RMS value of Centroid parameters
        end
        function [obj] = calcDelta(obj)
            x = abs(obj.fitParameters(:,2)-obj.Mean(:,2)); % absolute value from the mean in x
            y = abs(obj.fitParameters(:,4)-obj.Mean(:,4)); % absolute value from the mean in y
            r = sqrt(x.^2+y.^2); % pythag to find radius
            obj.delta = r; % assign value to obj
        end
        %=============================================%
        % methods that plot data from object
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
            plot(obj.PSD(:,1),10*log10(obj.PSD(:,2)))
            xlabel('f (Hz)')
            ylabel('Power/ Frequency (dB/Hz)')
            set(gca,'FontSize',16)
            xlim([0.5 obj.frameRate/2])
        end
        %=============================================%
        % methods that plot data from object or find frames which can be plotted%
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
        %=============================================%
        % methods to save object%
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
        %=============================================%
        % functions that fit image data %
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
            
            if centroid(1)<10 || centroid(1)> s(1)-10 ||centroid(2)<10 || centroid(2)>s(2)-10
                initialValues = [0,0,0,0,0,0,0];
                flag = 1;
                
            else
                
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
                flag = 0;
            end
            
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
                lb = [0,0,0.1,0,0.1,0,0];
                ub = [x0(1)*1.2,inf,inf,inf,inf,0,inf];
                % call LSQ curve fit without theta
                [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunction,x0,xdata,frame,lb,ub,options);
                
            end
        end
        %=============================================%
        % functions that condition frame for reading in and trimming%
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
        %=============================================%
        % interactive filter function%
        function [ry] = ifilter(ix,iy,icenter,iwidth,ishape,imode,ifilt)
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
function g = ngaussian(x,pos,wid,n)
%  ngaussian(x,pos,wid) = peak centered on x=pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  Shape is Gaussian when n=1, becomes more rectangular as n increases.
% Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
g = exp(-((x-pos)./(0.6005615.*wid)) .^(2*round(n)));
end
