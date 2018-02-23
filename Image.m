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
            obj.Mean = mean(obj.fitParameters); % Mean value of Centroid parameters
        end
        function [obj] = calcRange (obj)
            obj.Range = mean(obj.fitParameters); % Mean value of Centroid parameters
        end
        function [obj] = calcRMS (obj)
            obj.RMS = std(obj.fitParameters); % RMS value of Centroid parameters
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
            plot(obj.PSD(:,1),obj.PSD(:,2).^2)
            xlabel('f (Hz)')
            ylabel('|P1(f)|^2/f (pix^2/Hz)')
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

