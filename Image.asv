classdef Image
    properties
        Frame
        ApproxCentroid
        Centroid
        Params
        Flag
        Info
        Time
        Frequency
        Scatter
        FFT
        Residuals
        Mean
        RMS
        Range
        FitResidual
        Ci
        Filename
        NoFrames
        FrameNumbers
        PeakFrq
        cuts
        step
        Rho
        Name
        PixelSize
        stored_frames
    end
    methods
        function [obj] = RoughPeak(obj)
            for ii = 1:size(obj.Frame,3)
                frame=obj.Frame(:,:,ii);
                
                [peakamp,x_coord,y_coord] = nearest_pixel_centroid (frame);
                
                obj.ApproxCentroid(ii,1) = x_coord;
                obj.ApproxCentroid(ii,2) = y_coord;
                obj.ApproxCentroid(ii,3) = peakamp;
            end
        end
        function [obj] = SubNoise(obj)
            for ii = 1:size(obj.Frame,3)
                data=obj.Frame(:,:,ii);
                x_coord = obj.ApproxCentroid(1,1,ii);
                y_coord = obj.ApproxCentroid(1,2,ii);
                
                %%%%%%%%%%%----------------------noise reduction----------------------%%%%%%%%%%%
                
                samp = 0.1;
                
                %split the frame into 4 quadrants, br = bottom right, tr = top right,
                %bl = bottom left, tl = top left.
                bfindbr= [floor(size(data,2)*(1-samp)), size(data,2),floor(size(data,1)*(1-samp)), size(data,1)];
                bfindtr= [floor(size(data,2)*(1-samp)), size(data,2),1,floor(samp*size(data,1))];
                bfindbl= [1, floor(samp*size(data,2)),floor(size(data,1)*(1-samp)), size(data,1)] ;
                bfindtl= [1, floor(samp*size(data,2)),1,floor(samp*size(data,1))] ;
                
                % calculate the mean values in each quadrant
                br = mean(mean(data(bfindbr(3):bfindbr(4),bfindbr(1):bfindbr(2))));
                tr = mean(mean(data(bfindtr(3):bfindtr(4),bfindtr(1):bfindtr(2))));
                bl = mean(mean(data(bfindbl(3):bfindbl(4),bfindbl(1):bfindbl(2))));
                tl = mean(mean(data(bfindtl(3):bfindtl(4),bfindtl(1):bfindtl(2))));
                
                % Find the location of the PSF and use other 3 quadrants for background
                % noise reduction
                if y_coord > (size(data,1)/2)
                    if x_coord > (size(data,2)/2)
                        %disp('botright')
                        noise= mean([tr,bl,tl]);
                    else
                        %disp('botleft')
                        noise=mean([tr,br,tl]);
                    end
                    
                else
                    if x_coord>(size(data,2)/2);
                        %disp('topleft')
                        noise= mean([tr,br,bl]);
                    else
                        %disp('topright')
                        noise= mean([tl,br,bl]);
                    end
                end
                
                
                data = data-noise; % subtract noise of all data
                
                obj.Frame(:,:,ii) = data-noise;
            end
        end
        function [obj] = FinePeak(obj,FitForOrientation)
            
            % How many pixels do you want the Gaussian fit to extend over? (i.e. how big
            % is your PSF?) The framesize heavily impacts Gaussian width parameters.
            % Think about this one the most before moving on. Different framesizes will
            % lead to different function parameters depending on the noise, aberrations,
            % speckles outside the core region if they are included.
            
            for ii = 1:size(obj.Frame,3)
                
                x_coord = obj.ApproxCentroid(ii,1);
                y_coord = obj.ApproxCentroid(ii,2);
                peakamp = obj.ApproxCentroid(ii,3);
                frame = obj.Frame(:,:,ii);
                
                %x0 = [Amp,xo,wx,yo,wy,fi];
                x0 = [peakamp,x_coord,10,y_coord,10,0];
                
                %Fit the peak
                [x,flag] = sub_pixel_peak(frame,x0,FitForOrientation);
                
                % Assign object properties
                obj.Centroid(ii,1) = x(2); %xcen
                obj.Centroid(ii,2) = x(4); %ycen
                obj.Centroid(ii,3) = x(1); %peak
                obj.Centroid(ii,4) = x(3); %x width
                obj.Centroid(ii,5) = x(5); %y width
                obj.Centroid(ii,6) = x(6); %rotation
                obj.Flag = flag;
            end
        end
        function [obj] = Process_to_object(obj,data)
            
            %%Data Format:
            %Column 1: Frame time microseconds
            %Column 2: Actual time (based upon first frame timestamp to nearest second)
            %Column 3: X sub-pixel peak (Gaussian fit)
            %Column 4: Y sub-pixel peak (Gaussian fit)
            %Column 5: X sigma
            %Column 6: Y sigma
            %Column 7: theta, rotation angle of fit
            %Column 8: fit amplitude
            %Cloumn 9: constant offset fit (backgound)
            %Column 10: x centoid guess (variable used in script is x_coord)
            %Column 11: y centroid guess (variable used in script is y_coord)
            %Column 12: inital guess amp
            %Column 13: initial guess theta
            %Column 14: x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
            %Column 15: y_sigma
            %Column 16: Bad frame
            %Column 17: Frame Numbers
            %Column 18: FitResiduals
            
            obj.Time = data(:,1);
            obj.Centroid(:,1) = data(:,3);
            obj.Centroid(:,2) = data(:,4);
            obj.Centroid(:,3) = data(:,8);
            obj.Centroid(:,4) = data(:,5);
            obj.Centroid(:,5) = data(:,6);
            obj.Centroid(:,6) = data(:,7);
            obj.Centroid(:,7) = data(:,9);
            
            obj.ApproxCentroid(:,1) = data(:,10);
            obj.ApproxCentroid(:,2) = data(:,11);
            obj.FrameNumbers = data(:,17);
            obj.Flag = data(:,16);
            obj.FitResidual = data(:,18); %Fit residudals
            
        end
        function [obj] = Fourier_Transform (obj,xdata)
            X = xdata;
            Y = fft(X);
            L = length(X);
            Fs = obj.Frequency;
            % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);
            
            N = length(X);
            xdft = fft(X);
            xdft = xdft(1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            
            % Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.
            f = Fs*(0:(L/2))/L;
            obj.FFT(:,1) = f;
            obj.FFT(:,2) = P1;
        end
        function [obj,A] = Basler_Frames(obj,dataLocation,FitForOrientation)
            
            % N is the number of directories we want to process i.e Folders, Subfolders
            % M is the number of files in each directory
            
            % pathnames = strcat(dataStorageLocation,folderNames);
            % N = size(folderNames,2);
            
            file_id = strcat(dataLocation,'\*.tiff');
            allFiles = dir(file_id);
            allNames = {allFiles.name};
            M = size(allNames,2); %M is the number of files in the directory
            A = []; %reserved for frame data
            background = double((imread ('P:\iLocater\NA Spot Size test\Focus\Pair\Andrew\Background\acA3800-14um__21725217__20171101_144453723_0001.tiff')));
            for jj= 1:M
                allLocations{jj,1} = strcat(dataLocation,'\',allNames(jj));
                data = double(imread(char(allLocations{jj,1})));
                %                 data = data-background;
                
                %                 x0 = [peakamp,x_coord,10,y_coord,10,0]; %x0 = [Amp,xo,wx,yo,wy,fi];
                %                 [x,flag,x0,rms] = sub_pixel_peak_V2(data(:,:,jj),x0,FitForOrientation);
                sigma = [4,4];
                
                %                 if jj == 1
                [amp,x_coord,y_coord] = nearest_pixel_centroid (data(:,:));
                delta = 10;
                if jj == 1
                    left = max(floor(x_coord-delta*sigma(1)),1);
                    right = min(floor(x_coord+delta*sigma(1)),length(data(1,:)));
                    top = max(floor(y_coord-delta*sigma(2)),1);
                    bot = min(floor(y_coord+delta*sigma(2)),length(data(:,1)));
                    obj.cuts = [left,right,top,bot];
                end
                x0 = [amp,x_coord,sigma(1),y_coord,sigma(2),0,0]; %x0 = [Amp,xo,wx,yo,wy,theta,offset];
                %                 end
                
                %fit to obtain locations (especially good for largeor oversampled images)
                %                 [sigma,amp,cuts(:,1),center] = coarsefit(data);
                
                %2D fitting using least squares
                [x,flag,x0,residual,fiterr] = sub_pixel_peak_V2(data,x0,FitForOrientation,obj.cuts);
                
                if max(max(data)) > 254 || x(1) > 254%saturation limit
                    flag = 1;
                elseif x(2) < 0.1 || x(3) < 0.1% psf way too small
                    flag = 1;
                end
                
                frame(1) = NaN;
                frame(2) = NaN;
                frame(3) = x(2); % x(2) = x_centroid
                frame(4) = x(4); % x(4) = y_centroid
                frame(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
                frame(6) = x(5);  % x(5) = y_sigma
                frame(7) = x(6);  % x(6) = theta, rotation angle in radians
                frame(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
                frame(9) = x(7); % constant offset fit (i.e backgound)
                
                frame(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
                frame(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
                frame(12) = x0(1); % inital guess amp
                frame(13) = x0(3); % inital guess x wdith
                frame(14) = x0(5); % initial guess y width
                frame(15) = x0(6); % initial guess theta
                
                frame(16) = flag; % flag for bad frames
                frame(17) = 1; %FRAME NUMBER
                frame(18) = residual; %Fit residual
                obj.Ci(:,jj) = fiterr; %Fit error
                obj.Frame(:,:,jj) = data;
                A = [A; frame];
            end
            
            %             obj.Frame = data;
        end
        function [obj,frame] = SBIG_Frames(obj,filename,FitForOrientation)
            
            [data]=fitsread(filename);
            file = 'Background';
            ext = '.fit';
            path = 'P:\iLocater\GSFC_Fibers\NA2\1050nm\9_1\';
            darkframe = [path,file,ext];
            
            [background] = fitsread(darkframe);
            %[peakamp,x_coord,y_coord] = nearest_pixel_centroid (data);
            data = data-background;
            data =data./max(max(data));
            %fit to obtain locations (especially good for large images
            [sigma,amp,cuts(:,1),center] = coarsefit(data);
            %initial guess values
            x0 = [amp,center(2),sigma(2),center(1),sigma(1),0,0]; %x0 = [Amp,xo,wx,yo,wy,theta,offset];
            %2D fitting using least squares
            [x,flag,x0,residual,fiterr] = sub_pixel_peak_V2(data,x0,FitForOrientation,cuts);
            
            obj.cuts = cuts;
            
            frame(1) = NaN;
            frame(2) = NaN;
            frame(3) = x(2); % x(2) = x_centroid
            frame(4) = x(4); % x(4) = y_centroid
            frame(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
            frame(6) = x(5);  % x(5) = y_sigma
            frame(7) = x(6);  % x(6) = theta, rotation angle in radians
            frame(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
            frame(9) = x(7); % constant offset fit (i.e backgound)
            
            frame(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
            frame(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
            frame(12) = x0(1); % inital guess amp
            frame(13) = x0(3); % inital guess x wdith
            frame(14) = x0(5); % initial guess y width
            frame(15) = x0(6); % initial guess theta
            
            frame(16) = flag; % flag for bad frames
            frame(17) = 1; %FRAME NUMBER
            frame(18) = residual; %Fit residual
            obj.Ci = fiterr; %Fit error
            obj.Frame = data;
            
        end
        function [obj,frame] = Sim_Frames(obj,filename,FitForOrientation)
            
            load(filename)
            [data]=(PSF(:,:,3));
            obj.PixelSize = [mean(diff(PSF(1,:,1))),mean(diff(PSF(:,1,2)))];
            %[peakamp,x_coord,y_coord] = nearest_pixel_centroid (data);
            data =data./max(max(data));
            %fit to obtain locations (especially good for large images
            [sigma,amp,Cuts(:,1),center] = coarsefit(data);
            %initial guess values
            x0 = [amp,center(2),sigma(2),center(1),sigma(1),0,0]; %x0 = [Amp,xo,wx,yo,wy,theta,offset];
            %2D fitting using least squares
            [x,flag,x0,residual,fiterr] = sub_pixel_peak_V2(data,x0,FitForOrientation,Cuts);
            
            obj.cuts = Cuts;
            
            frame(1) = NaN;
            frame(2) = NaN;
            frame(3) = x(2); % x(2) = x_centroid
            frame(4) = x(4); % x(4) = y_centroid
            frame(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
            frame(6) = x(5);  % x(5) = y_sigma
            frame(7) = x(6);  % x(6) = theta, rotation angle in radians
            frame(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
            frame(9) = x(7); % constant offset fit (i.e backgound)
            
            frame(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
            frame(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
            frame(12) = x0(1); % inital guess amp
            frame(13) = x0(3); % inital guess x wdith
            frame(14) = x0(5); % initial guess y width
            frame(15) = x0(6); % initial guess theta
            
            frame(16) = flag; % flag for bad frames
            frame(17) = 1; %FRAME NUMBER
            frame(18) = residual; %Fit residual
            obj.Ci = fiterr; %Fit error
            obj.Frame = data;
        end
        function [obj] = SMF(obj)
            [X,Y] = meshgrid(1:size(obj.Frame,2),1:size(obj.Frame,1));
            %             center = [obj.Mean(1), obj.Mean(2)];
            Sigma = [3.6,3.6];
            for ii = 1:size(obj.Frame,3)
                center = [obj.Centroid(ii,1), obj.Centroid(ii,2)];
                Fiber = circ_gauss(X,Y,Sigma,center);
                Product = Fiber.*conj(obj.Frame(:,:,ii));
                D1 = abs(obj.Frame(:,:,ii)).^2;
                D2 = abs(Fiber).^2;
                num = abs(sum(sum(Product)).^2);
                den = abs(sum(sum(D1))*sum(sum(D2)));
                
                %Outputs
                obj.Rho(ii) = num/den;
            end
            
        end
        function [] = PeakPlot(obj,data_number)
            InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
            if nargin <2
                data_number = 1;
            end
            
            ii= obj.stored_frames(data_number);
            
            data = obj.Frame(:,:,data_number);
            x(2) = obj.Centroid(ii,1); %xcen
            x(4) = obj.Centroid(ii,2); %ycen
            x(1) = obj.Centroid(ii,3); %peak
            x(3) = obj.Centroid(ii,4); %x width
            x(5) = obj.Centroid(ii,5); %y width
            x(6) = obj.Centroid(ii,6); %rotation
            x(7) = obj.Centroid(ii,7); %constant
            dispangle = x(6)*180/pi;
            
            if isempty(obj.cuts) == 1
                [X,Y] = meshgrid(1:size(data,2),1:size(data,1));
                Z = data;
            else
                cuts = obj.cuts;
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
            
            
            %             string1 = ['Amplitude','    X-Center', '    X-Width','    Y-Center',...
            %                 '    Y-Width','    Angle','    Constant'];
            %             string3 = [num2str(x(1),'%100.3f'),'   ',num2str(x(2),'%100.3f'),'    ',num2str(x(3),'% 100.3f')...
            %                 ,'     ',num2str(x(4),'%100.3f'),'    ',num2str(x(5),'%100.3f'),'        ',num2str(dispangle,'%100.3f')...
            %                 ,'    ',num2str(x(7),'%100.3f')];
            %             text(min(X(1,:)),min(Y(:,1))-0.12*(max(Y(:,1))-min(Y(:,1))),string1,'Color','k','fontsize',12)
            %             text(min(X(1,:)),min(Y(:,1))-0.12*(max(Y(:,1))-min(Y(:,1)))-1.5,string3,'Color','red','fontsize',12)
            %
            
            
        end
        function [] = HistPlot(obj)
            
            %figure creation
            hf2 = figure();
            set(hf2, 'Position', [20 20 950 950])
            alpha(0)
            p1 = subplot(4,4,[5,6,7,9,10,11,13,14,15]);
            plot(obj.Centroid(:,1),obj.Centroid(:,2),'.')
            xlim([obj.Mean(1,1)-4*obj.RMS(1,1), obj.Mean(1,1)+4*obj.RMS(1,1)])
            ylim([obj.Mean(1,2)-4*obj.RMS(1,2), obj.Mean(1,2)+4*obj.RMS(1,2)])
            xlabel('Horizontal Pixels')
            ylabel('Vertical Pixels')
            axis equal
            box on
            p1.FontSize = 12;
            hold on
            
            %Histogram plots
            p2 = subplot(4,4,(1:3));
            histogram(obj.Centroid(:,1),round(sqrt(length(obj.Centroid(:,1)))))
            xlim([obj.Mean(1,1)-4*obj.RMS(1,1), obj.Mean(1,1)+4*obj.RMS(1,1)])
            p2.FontSize = 12;
            
            p3 = subplot(4,4,[8,12,16]);
            histogram(gca,obj.Centroid(:,2),round(sqrt(length(obj.Centroid(:,1)))))
            xlim([obj.Mean(1,2)-4*obj.RMS(1,2), obj.Mean(1,2)+4*obj.RMS(1,2)])
            p3.FontSize = 12;
            set(gca,'view',[90 -90])
            figure(gcf) % bring current figure to front
            hold off
            
            subplot(4,4,4)
            text(0,+0.5,{['Scatter Statistics:'],['STD (\sigma_x): ',num2str(obj.RMS(1),'%100.2f'),' pix'],['STD (\sigma_y): ',num2str(obj.RMS(2),'%100.2f'),' pix'],...
                ['Range x: ',num2str(obj.Range(1),'%100.2f'),' pix'],['Range y: ',num2str(obj.Range(2),'%100.2f'),' pix']},'FontSize',14);
            axis off
            
        end
        function [] = AngHistPlot(obj)
            
            %fiber
            pixscl=2.33; %mas/pixel
            %pixscl=0.604*6.5; % mas/um *6.4um/pix = mas/pix NIC focal plane with ANDOR
            fiber = 1/pixscl; %40/1.7/pixscl;
            dispx = (obj.Mean(1,1)-obj.Centroid(:,1));
            dispy = (obj.Mean(1,2)-obj.Centroid(:,2));
            
            %figure creation
            hf2 = figure();
            set(hf2, 'Position', [20 20 950 950])
            alpha(0)
            subplot(4,4,[5,6,7,9,10,11,13,14,15])
            plot(dispx./fiber,dispy./fiber,'.')
            xlim([-4*obj.RMS(1,1)/fiber, 4*obj.RMS(1,1)/fiber])
            %             ylim([-4*obj.RMS(1,2)/fiber, 4*obj.RMS(1,2)/fiber])
            ylim(xlim)
            xlabel('\Delta x-centroid (mas)')
            ylabel('\Delta y-centroid (mas)')
            set(gca,'FontSize',16)
            axis equal
            box on
            hold on
            
            %Histogram plots
            subplot(4,4,(1:3))
            histogram(dispx/fiber,round(sqrt(length(obj.Centroid(:,1)))))
            xlim([-4*obj.RMS(1,1)/fiber, 4*obj.RMS(1,1)/fiber])
            set(gca,'FontSize',18)
            
            subplot(4,4,[8,12,16])
            histogram(gca,dispy/fiber,round(sqrt(length(obj.Centroid(:,1)))))
            xlim([-4*obj.RMS(1,2)/fiber, 4*obj.RMS(1,2)/fiber])
            set(gca,'FontSize',16)
            set(gca,'view',[90 -90])
            figure(gcf) % bring current figure to front
            hold off
            
            
            %Statistics
            subplot(4,4,4)
            text(0,+0.5,{['Scatter Statistics:'],['STD (\sigma_x): ',num2str(obj.RMS(1)*pixscl,'%100.2f'),' mas'],['STD (\sigma_y): ',num2str(obj.RMS(2)*pixscl,'%100.2f'),' mas'],...
                ['Range x: ',num2str(obj.Range(1)*pixscl,'%100.2f'),' mas'],['Range y: ',num2str(obj.Range(2)*pixscl,'%100.2f'),' mas']},'FontSize',14);
            axis off
            
        end
        function [] = FFTPlot(obj)
            hf = figure();
            set(hf, 'Position', [20 20 1300 350])
            plot(obj.FFT(:,1),obj.FFT(:,2).^2)
            %             title('Single-Sided Amplitude Spectrum of X(t)')
            xlabel('F (Hz)')
            ylabel('|P1(f)|')
            set(gca,'FontSize',16)
            xlim([0.5 obj.Frequency/2])
        end
    end
    methods(Static)
        function [frame_time,image_data,cuts] = img_process(andor_file,lowMem,FitForOrientation,frames,mem_step)
            
            image_data = [];
            %%Data Format:
            %Column 1: Frame time microseconds
            %Column 2: Actual time (based upon first frame timestamp to nearest second)
            %Column 3: X sub-pixel peak (Gaussian fit)
            %Column 4: Y sub-pixel peak (Gaussian fit)
            %Column 5: X sigma
            %Column 6: Y sigma
            %Column 7: theta, rotation angle of fit
            %Column 8: fit amplitude
            %Cloumn 9: constant offset fit (backgound)
            %Column 10: x centoid guess (variable used in script is x_coord)
            %Column 11: y centroid guess (variable used in script is y_coord)
            %Column 12: inital guess amp
            %Column 13: initial guess theta
            %Column 14: x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
            %Column 15: y_sigma
            %Column 16: Bad frame
            
            frame_limit = 0; %If 0 or less, reads all frames. Otherwise only reads this many frames.
            %exclude_list = [31930,33201,37237,50001]; %NuCom_sub400
            %exclude_list = [19662,19721:19733,100001]; %Australis_6
            exclude_list = [];
            
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
                        
                        A =[];
                        for ii = (frames)
                            % Below grabs a Frame
                            [rc,data]=atsif_getframe(signal,ii-1,and_size);
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
                                
                                %This shapes a frame
                                width = ((right - left)+1)/hBin;
                                height = ((top-bottom)+1)/vBin;
                                newdata=double(reshape(data,width,height));
                                newdata = flipud(rot90(newdata));
                                if ii == 1
                                    image_data(:,:,1) = newdata;
                                end
                                
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
                                frame_time(2) = time_zero + frame_time(1) / (24 * 60 * 60 * 1e6) - (10/24);% calculate the matlab time
                                
                                %if(isempty(find(exclude_list == ii)))
                                %if lowMem == 1
                                %[peakamp,x_coord,y_coord] = nearest_pixel_centroid (newdata);
                                %x0 = [peakamp,x_coord,10,y_coord,10,0]; %x0 = [Amp,xo,wx,yo,wy,fi];
                                %[x,flag,x0,residual] = sub_pixel_peak(newdata,x0,FitForOrientation);
                                
                                %newdata = newdata - newdata(1,1);
                                %                                 if lowMem ==1
                                if ii == 1
                                    [sigma,amp,cuts(:,ii),center] = coarsefit(newdata);
                                else
                                    [amp,x_coord,y_coord] = nearest_pixel_centroid (newdata);
                                    center = [y_coord,x_coord];
                                    cuts(:,ii) = cuts(:,1);
                                    sigma(1) = A(1,6);
                                    sigma(2) = A(1,5);
                                end
                                %                                 else
                                %                                     [sigma,amp,cuts(:,ii),center] = coarsefit(newdata);
                                %                                 end
                                
                                x0 = [amp,center(2),sigma(2),center(1),sigma(1),0,0]; %x0 = [Amp,xo,wx,yo,wy,theta,offset];
                                [x,flag,x0,residual] = sub_pixel_peak_V2(newdata,x0,FitForOrientation,cuts(:,ii));
                                
                                frame_time(3) = x(2); % x(2) = x_centroid
                                frame_time(4) = x(4); % x(4) = y_centroid
                                frame_time(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
                                frame_time(6) = x(5);  % x(5) = y_sigma
                                frame_time(7) = x(6);  % x(6) = theta, rotation angle in radians
                                frame_time(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
                                frame_time(9) = x(7); % constant offset fit (i.e backgound)
                                
                                frame_time(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
                                frame_time(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
                                frame_time(12) = x0(1); % inital guess amp
                                frame_time(13) = x0(3); % inital guess x wdith
                                frame_time(14) = x0(5); % initial guess y width
                                frame_time(15) = x0(6); % initial guess theta
                                
                                frame_time(16) = flag; % flag for bad frames
                                frame_time(17) = ii-1; %FRAME NUMBER
                                frame_time(18) = residual; %Fit residual
                                
                                %                                     end
                                
                                %                                 else
                                %                                     frame_time(3:15) = NaN;
                                %                                     frame_time(16) = 1;
                                %                                 end
                                
                                if(mod(ii,50) == 0)
                                    fprintf('%i\n',ii)
                                end
                                
                            end
                            A = [A; frame_time];
                            if ii>1
                                if lowMem == 0 % frames are stored in image_data if there is memory
                                    jj = 1;
                                    image_data = cat(3,image_data,newdata);
                                    stored_frames(jj) = 1;
                                else
                                    if rem(ii,mem_step) == 0
                                        jj = jj+1;
                                        image_data = cat(3,image_data,newdata);
                                        stored_frames(jj) = ii;
                                    end
                                    
                                end
                            end
                            
                        end
                        
                        frame_time = A; %reassign dummy variable to frame time again
                    end
                end
                atsif_closefile;
            else
                disp('Could not load file.  ERROR - ');
                disp(rc);
            end
            
            %Calculating displacements and average pixel location
            %             row_count = 0;
            %             if lowMem == 1;
            %             %Remove bad images first so not to screw things up
            %             row_remove = find(isnan(frame_time(:,6)));
            %             row_count = row_count + size(row_remove);
            %             frame_time = removerows(frame_time,row_remove);
            %             row_remove = find(isnan(frame_time(:,7)));
            %             row_count = row_count + size(row_remove);
            %             frame_time = removerows(frame_time,row_remove);
            %             end
            
            %             frame_time = frame_time(1:10,:);
        end
        
        
        function [frame_time,image_data,cuts,stored_frames] = img_process_V2(andor_file,lowMem,FitForOrientation,frames,mem_step)
            
            image_data = [];
            %%Data Format:
            %Column 1: Frame time microseconds
            %Column 2: Actual time (based upon first frame timestamp to nearest second)
            %Column 3: X sub-pixel peak (Gaussian fit)
            %Column 4: Y sub-pixel peak (Gaussian fit)
            %Column 5: X sigma
            %Column 6: Y sigma
            %Column 7: theta, rotation angle of fit
            %Column 8: fit amplitude
            %Cloumn 9: constant offset fit (backgound)
            %Column 10: x centoid guess (variable used in script is x_coord)
            %Column 11: y centroid guess (variable used in script is y_coord)
            %Column 12: inital guess amp
            %Column 13: initial guess theta
            %Column 14: x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
            %Column 15: y_sigma
            %Column 16: Bad frame
            
            frame_limit = 0; %If 0 or less, reads all frames. Otherwise only reads this many frames.
            %exclude_list = [31930,33201,37237,50001]; %NuCom_sub400
            %exclude_list = [19662,19721:19733,100001]; %Australis_6
            exclude_list = [];
            
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
                        
                        A =[];
                        for ii = (frames)
                            % Below grabs a Frame
                            [rc,data]=atsif_getframe(signal,ii-1,and_size);
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
                                
                                %This shapes a frame
                                width = ((right - left)+1)/hBin;
                                height = ((top-bottom)+1)/vBin;
                                newdata=double(reshape(data,width,height));
                                newdata = flipud(rot90(newdata));
                                if ii == 1
                                    image_data(:,:,1) = newdata;
                                end
                                
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
                                frame_time(2) = time_zero + frame_time(1) / (24 * 60 * 60 * 1e6) - (10/24);% calculate the matlab time
                                
                                %if(isempty(find(exclude_list == ii)))
                                %if lowMem == 1
                                %[peakamp,x_coord,y_coord] = nearest_pixel_centroid (newdata);
                                %x0 = [peakamp,x_coord,10,y_coord,10,0]; %x0 = [Amp,xo,wx,yo,wy,fi];
                                %[x,flag,x0,residual] = sub_pixel_peak(newdata,x0,FitForOrientation);
                                
                                %newdata = newdata - newdata(1,1);
                                %                                 if lowMem ==1
                                if ii == 1
                                    jj = 1;
                                    [sigma,amp,cuts(:,ii),center] = coarsefit(newdata);
                                    stored_frames(1) = 1;
                                else
                                    [amp,x_coord,y_coord] = nearest_pixel_centroid (newdata);
                                    center = [y_coord,x_coord];
                                    cuts(:,ii) = cuts(:,1);
                                    sigma(1) = A(1,6);
                                    sigma(2) = A(1,5);
                                end
                                %                                 else
                                %                                     [sigma,amp,cuts(:,ii),center] = coarsefit(newdata);
                                %                                 end
                                
                                x0 = [amp,center(2),sigma(2),center(1),sigma(1),0,0]; %x0 = [Amp,xo,wx,yo,wy,theta,offset];
                                [x,flag,x0,residual] = sub_pixel_peak_V2(newdata,x0,FitForOrientation,cuts(:,ii));
                                
                                frame_time(3) = x(2); % x(2) = x_centroid
                                frame_time(4) = x(4); % x(4) = y_centroid
                                frame_time(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
                                frame_time(6) = x(5);  % x(5) = y_sigma
                                frame_time(7) = x(6);  % x(6) = theta, rotation angle in radians
                                frame_time(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
                                frame_time(9) = x(7); % constant offset fit (i.e backgound)
                                
                                frame_time(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
                                frame_time(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
                                frame_time(12) = x0(1); % inital guess amp
                                frame_time(13) = x0(3); % inital guess x wdith
                                frame_time(14) = x0(5); % initial guess y width
                                frame_time(15) = x0(6); % initial guess theta
                                
                                frame_time(16) = flag; % flag for bad frames
                                frame_time(17) = ii-1; %FRAME NUMBER
                                frame_time(18) = residual; %Fit residual
                                
                                %                                     end
                                
                                %                                 else
                                %                                     frame_time(3:15) = NaN;
                                %                                     frame_time(16) = 1;
                                %                                 end
                                
                                if(mod(ii,50) == 0)
                                    fprintf('%i\n',ii)
                                end
                                
                            end
                            A = [A; frame_time];
                            
                            if ii>1
                                if lowMem == 0 % frames are stored in image_data if there is memory
                                    
                                    image_data = cat(3,image_data,newdata);
                                    
                                else
                                    if rem(ii,mem_step) == 0
                                        jj = jj+1;
                                        image_data = cat(3,image_data,newdata);
                                        stored_frames(jj) = ii;
                                    end
                                    
                                end
                            end
                            
                        end
                        
                        frame_time = A; %reassign dummy variable to frame time again
                    end
                end
                atsif_closefile;
            else
                disp('Could not load file.  ERROR - ');
                disp(rc);
            end
            
            %Calculating displacements and average pixel location
            %             row_count = 0;
            %             if lowMem == 1;
            %             %Remove bad images first so not to screw things up
            %             row_remove = find(isnan(frame_time(:,6)));
            %             row_count = row_count + size(row_remove);
            %             frame_time = removerows(frame_time,row_remove);
            %             row_remove = find(isnan(frame_time(:,7)));
            %             row_count = row_count + size(row_remove);
            %             frame_time = removerows(frame_time,row_remove);
            %             end
            
            %             frame_time = frame_time(1:10,:);
        end 
        function [no_frames] = calc_num_frames(filename)
            frame_limit = 0; %If 0 or less, reads all frames. Otherwise only reads this many frames.
            rc=atsif_setfileaccessmode(0);
            rc=atsif_readfromfile(filename);
            if (rc == 22002)
                signal=0;
                [rc,present]=atsif_isdatasourcepresent(signal);
                if present
                    [rc,no_frames]=atsif_getnumberframes(signal);
                    
                    if(frame_limit > 0)
                        no_frames = frame_limit;
                    end
                    
                end
            end
        end
    end
end

function [sigma_guess,amp_guess,cuts_guess,center_guess] = coarsefit(frame)

% frame = frame./mean(mean(frame));
% frame = frame.^2;

delta = 2;
rsum = sum(frame,1)';
csum = sum(frame,2);

geqn = 'a.*exp(-0.5*((x-b)/c).^2)+d';

[rmax,rloc] = max(rsum);
[cmax,cloc] = max(csum);

rstartpoints = [rmax rloc 10 0];
cstartpoints = [cmax cloc 10 0];
%
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

% [rfit,rgof] = fit((1:length(rsum))',rsum,geqn,'start',rstartpoints,optionsR);
% [cfit,cgof] = fit((1:length(csum))',csum,geqn,'start',cstartpoints,options);
[rfit,rgof] = fit((1:length(rsum))',rsum,geqn,optionsR);
[cfit,cgof] = fit((1:length(csum))',csum,geqn,optionsC);


% assert(rgof.adjrsquare>=0.0001,'bad row fit');
% assert(cgof.adjrsquare>=0.0001,'bad column fit');

c=coeffvalues(cfit);
r=coeffvalues(rfit);

rcenter = r(2);
ccenter = c(2);

sigma = max(c(3),r(3));
% sigma = 20;

left = max(floor(rcenter-delta*sigma),1);
right = min(floor(rcenter+delta*sigma),length(rsum));
top = max(floor(ccenter-delta*sigma),1);
bot = min(floor(ccenter+delta*sigma),length(csum));
cuts_guess = [left,right,top,bot];

% center_guess = [ccenter,rcenter];
center_guess = [cloc,rloc];
sigma_guess =([c(3),r(3)]);
amp_guess = frame(cloc,rloc);
end
function [x,flag,x0,nrms,fiterror] = sub_pixel_peak_V2(frame,x0,FitForOrientation,cuts)
options=optimset('Diagnostics','off','Display','none'); % options set for LSQ fit

%%%%%%%%%%%----------------------Frame dimensions for fitting----------------------%%%%%%%%%%%

%framesize is user defined size of the area which will be used for Gaussian
%fitting. i.e ony choose a size (in pixels) which is/should be Gaussian.
%framesize will always produce a square frame.

% [X,Y] = meshgrid(cuts(1):cuts(2),size(frame,1)+1-cuts(3):-1:size(frame,1)+1-cuts(4));
[X,Y] = meshgrid(cuts(1):cuts(2),cuts(3):cuts(4));

% Meshgrid steps over '0' which makes the frame 1 pixel larger than
% desired. Truncates to correct size. Offsets all frame values (i.e pixels)
% to be centered at location of PSF - represents the actual detector location

% X = X(1:yradius*2,1:xradius*2)+floor(center(1,2));
% Y = Y(yradius*2:-1:1,1:xradius*2)+floor(center(1,1));

% xdata = zeros(size(X,1),size(Y,2),2); %creates xdata - map of search pixel space.
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
    %                     obj.FineCentroid(1,1,ii) = NaN; %xcen
    %                     obj.FineCentroid(1,2,ii) = NaN; %ycen
    %                     obj.FineCentroid(1,3,ii) = NaN; %peak
    %                     obj.FineCentroid(1,4,ii) = NaN; %x width
    %                     obj.FineCentroid(1,5,ii) = NaN; %y width
    %                     obj.FineCentroid(1,6,ii) = NaN; %rotation
    
    x= [NaN,NaN,NaN,NaN,NaN,NaN];
    flag = 1;
else
    %Jonathan End Edit
    
    Z = frame(cuts(3):cuts(4),cuts(1):cuts(2)); % cut data at the locations corresponding to the chosen frame size.
    
    %%%%%%%%%%%----------------------LSQ fitting for 2D fits----------------------%%%%%%%%%%%
    
    if nargin == 4
    else
        FitForOrientation =0;
    end
    
    if FitForOrientation == 1
        lb = [0,0,0,0,0,-pi/4,0];
        ub = [inf,inf,inf,inf,inf,pi/4,inf];
        % call LSQ curve fit with theta
        [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub,options);
    else
        x0(6) = 0;
        lb = [0,0,0,0,0,0,0];
        ub = [inf,inf,inf,inf,inf,0,inf];
        % call LSQ curve fit without theta
        [x,resnorm,residual,exitflag,output,lambda,J1] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub,options);
        
    end
    
    J1 = full(J1);
    rms = sqrt(mean(mean(residual.^2)));%rms value for the 2D fit
    nrms = rms/(max(max(Z))-min(min((Z))));%normalizd rms value for the 2D fit
    parci = nlparci(x,residual,'Jacobian',J1);
    fiterror = (diff(parci,[],2)./2)';
    
    return
    
    %remap my parameters to work with this method (mostly need this for
    %Jacobian to work (The one from fit output gives NAN when used to find
    %Confidence intervals
    
    [xData, yData, zData] = prepareSurfaceData(X,Y,Z);
    xyData = {xData,yData};
    par(7) = x(7); % constant offset;
    par(1) = x(1); % amplitude
    par(5) = x(2); % x-center
    par(6) = x(3); % y-center
    par(3) = sqrt(2)*x(3); % x-width
    par(4) = sqrt(2)*x(5); % y-widtch
    par(2) = x(6); % theta
    [fiterr] = gaussian2Duncert(par,residual,xyData); %calculate confidence intervals
    
    %remap confidence itervals (fit parameters to match my ordering for
    %convenience
    fiterror(1) = fiterr(1); %amplitude
    fiterror(2) = fiterr(5); %x-center
    fiterror(3) = fiterr(6); %y-center
    fiterror(4) = fiterr(3)/sqrt(2); %x-sigma
    fiterror(5) = fiterr(4)/sqrt(2); %y-sigma
    fiterror(6) = fiterr(2); %theta
    fiterror(7) = fiterr(7); %constant
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
function [x_shift,y_shift] = xcorr_shift (frame,ref)
%%%%%%%%%%%%%%%%%%%%---------------Locate pixel shift-----------------%%%%%%%%%%%%%%%%%%%%%%%
C = normxcorr2(ref,frame);
m =size(ref,1);
n =size(ref,2);

[~,I] = max(frame(:)); %find max of linearized data
[I1,I2] = ind2sub(size(frame),I); % find x,y index locations of maximum count using linear index

x_shift = I1-(m-1)/2;
y_shift = I2-(m-1)/2;
end
function [peakamp,x_coord,y_coord] = nearest_pixel_centroid (frame)
%%%%%%%%%%%%%%%%%%%%---------------Locate pixel co-ordinates-----------------%%%%%%%%%%%%%%%%%%%%%%%
peakamp = max(max(frame)); % find maximum count pixel use amplitude as fit parameter
[~,I] = max(frame(:)); %find max of linearized data
[I1,I2] = ind2sub(size(frame),I); % find x,y index locations of maximum count using linear index
x_coord = I2; %guess for x_coord centroid
y_coord = I1; %guess for y_coord centroid
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


% Example:
%          F=circ_gauss(MatX,MatY,[1],[0 0]);
%          surface(F);
%--------------------------------------------------------------------------



SigmaX = Sigma(1);
SigmaY = Sigma(2);

X0 = center(1);
Y0 = center(2);

F = 1./(2.*pi.*SigmaX.*SigmaY) .*exp(-1./(2.).* ((X-X0).^2./SigmaX.^2 +(Y-Y0).^2./SigmaY.^2));



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
function [dpar] = gaussian2Duncert(par,resid,xy)
% get the confidence intervals
J = guassian2DJacobian(par,xy);
parci = nlparci(par,resid,'Jacobian',J);
dpar = (diff(parci,[],2)./2)';
%[zf,dzf] = nlpredci(@D2GaussFunction,xdata,par,resid,'Jacobian',J);
end
function J = guassian2DJacobian(par,xy)
% compute the jacobian
x = xy{1}; y = xy{2};
J(:,1) = exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2);
J(:,2) = -par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(3).^2 - (2.*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,3) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2)./par(3)^3;
J(:,4) = (2.*par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2)./par(4)^3;
J(:,5) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2 - (2.*sind(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2);
J(:,6) = par(1).*exp(- (cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))).^2./par(3).^2 - (cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))).^2./par(4).^2).*((2.*cosd(par(2)).*(cosd(par(2)).*(y - par(6)) - sind(par(2)).*(x - par(5))))./par(4).^2 + (2.*sind(par(2)).*(cosd(par(2)).*(x - par(5)) + sind(par(2)).*(y - par(6))))./par(3).^2);
J(:,7) = ones(size(x));
end
function rr = rsquared(z,zf,ze)
% reduced chi-squared
dz = z-zf;
rr = 1./(numel(z)-8).*sum(dz.^2./ze.^2); % minus 8 because there are 7 fit parameters +1 (DOF)
end
function z = gaussian2D(par,xy)
% compute 2D gaussian
z = par(7) + ...
    par(1)*exp(-(((xy{1}-par(5)).*cosd(par(2))+(xy{2}-par(6)).*sind(par(2)))./par(3)).^2-...
    ((-(xy{1}-par(5)).*sind(par(2))+(xy{2}-par(6)).*cosd(par(2)))./par(4)).^2);
end

%Unused/Outdated -- Do not move to next version!!

function [frame] = Basler_Frame_Grab(data)

N = size(data,2);

for ii = 1:N
    M = size(data(:,ii),1);
    
    for jj = 1:M
        frame=data{jj,ii};
        if isempty(frame)==0
        else
            disp('empty')
        end
        
    end
end
end
function [allLocations,data,A] = Dir_Calc(folderNames,dataLocation)

% N is the number of directories we want to process i.e Folders, Subfolders
% M is the number of files in each directory

% pathnames = strcat(dataStorageLocation,folderNames);
% N = size(folderNames,2);

for ii = 1:1 % we dont use the directory thing anymore.... 7/2017
    file_id = strcat(dataLocation,'\*.tiff');
    allFiles = dir(file_id);
    allNames = {allFiles.name};
    M = size(allNames,2); %M is the number of files in the directory
    A = [];
    for jj= 1:M
        allLocations{jj,ii} = strcat(pathnames{ii},'\',allNames(jj));
        data(:,:,jj) = double(imread(char(allLocations{jj,ii})));
        
        [peakamp,x_coord,y_coord] = nearest_pixel_centroid (data(:,:,jj));
        FitForOrientation = 1;
        x0 = [peakamp,x_coord,10,y_coord,10,0]; %x0 = [Amp,xo,wx,yo,wy,fi];
        [x,flag,x0] = sub_pixel_peak(data(:,:,jj),x0,FitForOrientation);
        
        frame(1) = NaN;
        frame(2) = NaN;
        frame(3) = x(2); % x(2) = x_centroid
        frame(4) = x(4); % x(4) = y_centroid
        frame(5) = x(3);  % x(3) = x_sigma, sigma is the width parameter (2.35*sigma = FHWM)
        frame(6) = x(5);  % x(5) = y_sigma
        frame(7) = x(6);  % x(6) = theta, rotation angle in radians
        frame(8) = x(1); % x(1) = max_count  peak count of the Gaussian fit
        frame(9) = x(7); % constant offset fit (i.e backgound)
        
        frame(10) = x0(2); % x0(2) = x centoid guess (variable used in script is x_coord)
        frame(11) = x0(4); % x0(3) = y centroid guess (variable used in script is y_coord)
        frame(12) = x0(1); % inital guess amp
        frame(13) = x0(3); % inital guess x wdith
        frame(14) = x0(5); % initial guess y width
        frame(15) = x0(6); % initial guess theta
        
        frame(16) = flag; % flag for bad frames
        frame(17) = ii-1; %FRAME NUMBER
        %         info = imfinfo(char(allLocations{jj,ii}));
        A = [A; frame];
    end
    
end

disp('finished!')

end
function [x1, y1] = ANDOR_find_peak(data)
%LMIRCam_PSF_Locate - Function to identify pixel position of peak in PSF.
%For LMIRCam, two values should be returned. If there is a significant
%discrepancy between the values of the PSFs, the one_psf_flag will be set.

%The algorithm tries to avoid hot pixels to find the approximate location
%of the brightest pixel.

peaks = 0;
frac_val = 0.30;
frac_change = 0.002;
frac_min = 0.002;

psf_flag = 0;

while (size(peaks,1) < 2) && (frac_val >= frac_min)
    peaks = FastPeakFind(data,max(max(data))*frac_val);
    frac_val = frac_val - frac_change;
end

x1 = peaks(1,1);
y1 = peaks(2,1);

end
function [cent, varargout]=FastPeakFind(d, thres, filt ,edg, res, fid)
% Analyze noisy 2D images and find peaks using local maxima (1 pixel
% resolution) or weighted centroids (sub-pixel resolution).
% The code is designed to be as fast as possible, so I kept it pretty basic.
% The code assumes that the peaks are relatively sparse, test whether there
% is too much pile up and set threshold or user defined filter accordingly.
%
% How the code works:
% In theory, each peak is a smooth point spread function (SPF), like a
% Gaussian of some size, etc. In reality, there is always noise, such as
%"salt and pepper" noise, which typically has a 1 pixel variation.
% Because the peak's PSF is assumed to be larger than 1 pixel, the "true"
% local maximum of that PSF can be obtained if we can get rid of these
% single pixel noise variations. There comes medfilt2, which is a 2D median
% filter that gets rid of "salt and pepper" noise. Next we "smooth" the
% image using conv2, so that with high probability there will be only one
% pixel in each peak that will correspond to the "true" PSF local maximum.
% The weighted centroid approach uses the same image processing, with the
% difference that it just calculated the weighted centroid of each
% connected object that was obtained following the image processing.  While
% this gives sub-pixel resolution, it can miss peaks that are very close to
% each other, and runs slightly slower. Read more about how to treat these
% cases in the relevant code commentes.
%
% Inputs:
% d     The 2D data raw image - assumes a Double\Single-precision
%       floating-point, uint8 or unit16 array. Please note that the code
%       casts the raw image to uint16 if needed.  If the image dynamic range
%       is between 0 and 1, I multiplied to fit uint16. This might not be
%       optimal for generic use, so modify according to your needs.
% thres A number between 0 and max(raw_image(:)) to remove  background
% filt  A filter matrix used to smooth the image. The filter size
%       should correspond the characteristic size of the peaks
% edg   A number>1 for skipping the first few and the last few 'edge' pixels
% res   A handle that switches between two peak finding methods:
%       1 - the local maxima method (default).
%       2 - the weighted centroid sub-pixel resolution method.
%       Note that the latter method takes ~20% more time on average.
% fid   In case the user would like to save the peak positions to a file,
%       the code assumes a "fid = fopen([filename], 'w+');" line in the
%       script that uses this function.
%
%Optional Outputs:
% cent        a 1xN vector of coordinates of peaks (x1,y1,x2,y2,...
% [cent cm]   in addition to cent, cm is a binary matrix  of size(d)
%             with 1's for peak positions. (not supported in the
%             the weighted centroid sub-pixel resolution method)
%
%Example:
%
%   p=FastPeakFind(image);
%   imagesc(image); hold on
%   plot(p(1:2:end),p(2:2:end),'r+')
%
%   Adi Natan (natan@stanford.edu)
%   Ver 1.7 , Date: Oct 10th 2013
%
%% defaults
if (nargin < 1)
    d=uint16(conv2(reshape(single( 2^14*(rand(1,1024*1024)>0.99995) ),[1024 1024]) ,fspecial('gaussian', 15,3),'same')+2^8*rand(1024));
    imagesc(d);
end

if ndims(d)>2 %I added this in case one uses imread (JPG\PNG\...).
    d=uint16(rgb2gray(d));
end

if isfloat(d) %For the case the input image is double, casting to uint16 keeps enough dynamic range while speeds up the code.
    if max(d(:))<=1
        d =  uint16( d.*2^16./(max(d(:))));
    else
        d = uint16(d);
    end
end

if (nargin < 2)
    thres = (max([min(max(d,[],1))  min(max(d,[],2))])) ;
end

if (nargin < 3)
    filt = (fspecial('gaussian', 7,3)); %if needed modify the filter according to the expected peaks sizes
end

if (nargin < 4)
    edg = 100;
end

if (nargin < 5)
    res = 1;
end

if (nargin < 6)
    savefileflag = false;
else
    savefileflag = true;
end

%% Analyze image
if any(d(:))  ; %for the case of non zero raw image
    
    d = medfilt2(d,[3,3]);
    
    % apply threshold
    if isa(d,'uint8')
        d=d.*uint8(d>thres);
    else
        d=d.*uint16(d>thres);
    end
    
    if any(d(:))   ; %for the case of the image is still non zero
        
        % smooth image
        d=conv2(single(d),filt,'same') ;
        
        % Apply again threshold (and change if needed according to SNR)
        d=d.*(d>0.9*thres);
        
        switch res % switch between local maxima and sub-pixel methods
            
            case 1 % peak find - using the local maxima approach - 1 pixel resolution
                
                % d will be noisy on the edges, and also local maxima looks
                % for nearest neighbors so edge must be at least 1. We'll skip 'edge' pixels.
                sd=size(d);
                [x y]=find(d(edg:sd(1)-edg,edg:sd(2)-edg));
                
                % initialize outputs
                cent=[];%
                cent_map=zeros(sd);
                
                x=x+edg-1;
                y=y+edg-1;
                for j=1:length(y)
                    if (d(x(j),y(j))>=d(x(j)-1,y(j)-1 )) &&...
                            (d(x(j),y(j))>d(x(j)-1,y(j))) &&...
                            (d(x(j),y(j))>=d(x(j)-1,y(j)+1)) &&...
                            (d(x(j),y(j))>d(x(j),y(j)-1)) && ...
                            (d(x(j),y(j))>d(x(j),y(j)+1)) && ...
                            (d(x(j),y(j))>=d(x(j)+1,y(j)-1)) && ...
                            (d(x(j),y(j))>d(x(j)+1,y(j))) && ...
                            (d(x(j),y(j))>=d(x(j)+1,y(j)+1));
                        
                        %All these alternatives were slower...
                        %if all(reshape( d(x(j),y(j))>=d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1))
                        %if  d(x(j),y(j)) == max(max(d((x(j)-1):(x(j)+1),(y(j)-1):(y(j)+1))))
                        %if  d(x(j),y(j))  == max(reshape(d(x(j),y(j))  >=  d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1))
                        
                        cent = [cent ;  y(j) ; x(j)];
                        cent_map(x(j),y(j))=cent_map(x(j),y(j))+1; % if a binary matrix output is desired
                        
                    end
                end
                
            case 2 % find weighted centroids of processed image,  sub-pixel resolution.
                % no edg requirement needed.
                
                % get peaks areas and centroids
                stats = regionprops(logical(d),d,'Area','WeightedCentroid');
                
                % find reliable peaks by considering only peaks with an area
                % below some limit. The weighted centroid method can be not
                % accurate if peaks are very close to one another, i.e., a
                % single peak will be detected, instead of the real number
                % of peaks. This will result in a much larger area for that
                % peak. At the moment, the code ignores that peak. If that
                % happens often consider a different threshold, or return to
                % the more robust "local maxima" method.
                % To set a proper limit, inspect your data with:
                % hist([stats.Area],min([stats.Area]):max([stats.Area]));
                % to see if the limit I used (mean+2 standard deviations)
                % is an appropriate limit for your data.
                
                rel_peaks_vec=[stats.Area]<=mean([stats.Area])+2*std([stats.Area]);
                cent=[stats(rel_peaks_vec).WeightedCentroid]';
                cent_map=[];
                
        end
        
        if savefileflag
            % previous version used dlmwrite, which can be slower than  fprinf
            %             dlmwrite([filename '.txt'],[cent],   '-append', ...
            %                 'roffset', 0,   'delimiter', '\t', 'newline', 'unix');+
            
            fprintf(fid, '%f ', cent(:));
            fprintf(fid, '\n');
            
        end
        
        
    else % in case image after threshold is all zeros
        cent=[];
        cent_map=zeros(size(d));
        if nargout>1 ;  varargout{1}=cent_map; end
        return
    end
    
else % in case raw image is all zeros (dead event)
    cent=[];
    cent_map=zeros(size(d));
    if nargout>1 ;  varargout{1}=cent_map; end
    return
end

%demo mode - no input to the function
if (nargin < 1); colormap(bone);hold on; plot(cent(1:2:end),cent(2:2:end),'rs');hold off; end

% return binary mask of centroid positions if asked for
if nargout>1 ;  varargout{1}=cent_map; end
end
function [x,flag,x0,nrms] = sub_pixel_peak(frame,x0,FitForOrientation)
framesize=200; %must be an EVEN number
options=optimset('Diagnostics','off','Display','none'); % options set for LSQ fit
% define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
yd = 500; %set how many pixels from the center guess the fit will seach in vertical direcion
xd = 500;
lsig = 1; %set minimum width parameter (1pixel is smallest)
usig = 1000; %set maximum width parameter (set larger than your PSF)
thetamin = 0; % sets the range of theta values. Should never need to change bounds on theta.
thetamax = pi/2;
xfactor = 2; %max fit value = amplitdue*xfactor, min fit amp = = amplitdue/xfactor
x_coord = x0(2);
y_coord = x0(4);
peakamp = x0(1);
x0(7) = 0;
%%%%%%%%%%%----------------------Frame dimensions for fitting----------------------%%%%%%%%%%%

%framesize is user defined size of the area which will be used for Gaussian
%fitting. i.e ony choose a size (in pixels) which is/should be Gaussian.
%framesize will always produce a square frame.


% Uses meshgrid to create a matrix with call values from -framesize to
% +framesize (see meshgrid funtion for more). This type of matrix is
% required for LSQcurvefit to define region of seach space

[X,Y] = meshgrid(-framesize/2:framesize/2,-framesize/2:framesize/2);

% Meshgrid steps over '0' which makes the frame 1 pixel larger than
% desired. Truncates to correct size. Offsets all frame values (i.e pixels)
% to be centered at location of PSF - represents the actual detector location
X = X(1:framesize,1:framesize)+x_coord;
Y = Y(1:framesize,1:framesize)+y_coord;

xdata = zeros(size(X,1),size(Y,2),2); %creates xdata - map of search pixel space.
xdata(:,:,1) = X; %layer 1 is X search space
xdata(:,:,2) = Y; %layer 2 is Y search space

Bot = Y(1,1);
Top = Y(framesize,1);
Left = X(1,1);
Right = X(1,framesize);

%%%%%%%%%%%----------------------remove bad frames----------------------%%%%%%%%%%%
flag = 0; % error check
if(Bot < 1 || Top > size(frame,1) || Left < 1 || Right > size(frame,2))
    %                     obj.FineCentroid(1,1,ii) = NaN; %xcen
    %                     obj.FineCentroid(1,2,ii) = NaN; %ycen
    %                     obj.FineCentroid(1,3,ii) = NaN; %peak
    %                     obj.FineCentroid(1,4,ii) = NaN; %x width
    %                     obj.FineCentroid(1,5,ii) = NaN; %y width
    %                     obj.FineCentroid(1,6,ii) = NaN; %rotation
    
    x= [NaN,NaN,NaN,NaN,NaN,NaN];
    flag = 1;
else
    %Jonathan End Edit
    
    Z = frame(Bot:Top,Left:Right); % cut data at the locations corresponding to the chosen frame size.
    
    %%%%%%%%%%%----------------------LSQ fitting----------------------%%%%%%%%%%%
    
    if nargin == 3
    else
        x0 = [peakamp,x_coord,9,y_coord,9,0]; %Inital guess parameters stored in array x0
        FitForOrientation =0;
    end
    
    if FitForOrientation == 1
        % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
        lb = [peakamp/xfactor,x_coord-xd,lsig,y_coord-yd,lsig,thetamin,0];
        ub = [peakamp*xfactor,x_coord+xd,usig,y_coord+yd,usig,thetamax,10000];
        % call LSQ curve fit with theta
        [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub,options);
        %         values = x;
    else
        x0(6) = 0;
        %         lb = [peakamp/xfactor,x_coord-xd,lsig,y_coord-yd,lsig,-inf,0];
        %         ub = [peakamp*xfactor,x_coord+xd,usig,y_coord+yd,usig,inf,10000];
        lb = [0,0,0,0,0,0,0];
        ub = [inf,inf,inf,inf,inf,0,inf];
        % call LSQ curve fit without theta
        [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub,options);
        %         values = x;
        
    end
    %         values(7)=x_coord; %assign to output varaible
    %         values(8)=y_coord;
    %         values(9)=residual;
    rms = sqrt(mean(mean(residual.^2)));
    nrms = rms/(max(max(Z))-min(min((Z))));
    
end
end
function [] = NormPeakPlot(obj)
for ii = 1:1;
    data = obj.Frame(:,:,ii);
    x_coord = obj.ApproxCentroid(1,1,ii);
    y_coord = obj.ApproxCentroid(1,2,ii);
    x(2) = obj.Centroid(1,1,ii); %xcen
    x(4) = obj.Centroid(1,2,ii); %ycen
    x(1) = obj.Centroid(1,3,ii); %peak
    x(3) = obj.Centroid(1,4,ii); %x width
    x(5) = obj.Centroid(1,5,ii); %y width
    x(6) = obj.Centroid(1,6,ii); %rotation
    dispangle = x(6)*180/pi;
    InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
    
    framesize = round(obj.Centroid(1,4,ii)*10);
    %                 if min(size(data))<100
    %                     edge=10*floor(min(size(obj.Frame))/10);
    %                     framesize = edge-min(floor(y_coord),floor(x_coord))-1;
    %                 else
    %                     framesize = 100;
    %                 end
    
    [X,Y] = meshgrid(-framesize/2:framesize/2);
    X = X (1:framesize,1:framesize);
    Y = Y (1:framesize,1:framesize);
    
    X = X+x_coord;
    Y = Y+y_coord;
    
    Bot = x_coord-framesize/2+1;
    Top = x_coord+framesize/2;
    Left = y_coord-framesize/2+1;
    Right = y_coord+framesize/2;
    Z = data(Left:Right,Bot:Top);
    
    
    %figure creation
    
    grid = 100;
    scale = 7;
    hf2 = figure(ii);
    set(hf2, 'Position', [0 0 scale*grid scale*grid])
    
    %                 axes same
    %                 axes equal
    
    subplot(10,10,[21:27,31:37,41:47,51:57,61:67,71:77,81:87])
    alpha(0)
    box on
    imagesc(X(1,:),Y(:,1)',Z)
    set(gca,'YDir','reverse')
    colormap('bone')
    %                 T = table(Age,Height,Weight,'RowNames',LastName);
    %                 uitable('Data',[obj.Params(ii,:);x],'RowName',{'Initial','Final'}, 'ColumnName', {'Peak', 'X', 'Y','Sigma X','Sigma Y','Theta'}, 'Position', scale*grid*[0.12 0.06 0.758 0.085]);
    %                                 string1 = ['       Amplitude','    X-Coordinate', '    X-Width','    Y-Coordinate','    Y-Width','     Angle'];
    %                 %             string2 = ['Set  ',num2str(xin(1),'% 100.3f'),'           ',num2str(xin(2),'% 100.3f'),'       ',num2str(xin(3),'% 100.3f'),'         ',num2str(xin(4), '% 100.3f'),'        ',num2str(xin(5), '% 100.3f'),'     ',num2str(dispangle, '% 100.3f')];
    %                 string3 = ['Fit   ',num2str(x(1),'% 100.3f'),'            ',num2str(x(2),'% 100.3f'),'      ',num2str(x(3),'% 100.3f'),'         ',num2str(x(4), '% 100.3f'),'        ',num2str(x(5), '% 100.3f'),'     ',num2str(dispangle, '% 100.3f')];
    %
    %                 text((x_coord)-framesize/3,+((y_coord)+framesize/2)+6,string1,'Color','k')
    %                 %             text((x_coord)-MdataSize/3,+((y_coord)+MdataSize/2)+8,string2,'Color','red')
    %                 text((x_coord)-framesize/3,+((y_coord)+framesize/2)+10,string3,'Color','red')
    
    % -----Calculate cross sections-------------
    %generate points along horizontal axis
    m = -tan(x(6));% Point slope formula
    b = (-m*x(2) + x(4));
    xvh = -framesize/2+x_coord:framesize/2+x_coord;
    yvh = xvh*m + b;
    hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
    %generate points along vertical axis
    mrot = -m;
    brot = (mrot*x(4) - x(2));
    yvv = -framesize/2+y_coord:framesize/2+y_coord;
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
    axis([-framesize/2+(x_coord) framesize/2+(x_coord) -framesize/2+(y_coord) framesize/2+(y_coord)])
    %
    
    ymin = 0;
    ymax = x(1);
    
    xdatafit = linspace(-framesize/2+(x_coord), framesize/2+(x_coord),300);
    hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
    
    ydatafit = linspace(-framesize/2+(y_coord), framesize/2+(y_coord),300);
    vdatafit = x(1)*exp(-(ydatafit-x(4)).^2/(2*x(5)^2));
    
    subplot(10,10,[1:7,11:17])
    xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
    plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
    axis([-framesize/2+(x_coord) framesize/2+(x_coord) ymin*1.1 ymax*1.5])
    
    subplot(10,10,[29:30,39:40,49:50,59:60,69:70,79:80,89:90])
    xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
    
    plot(vPoints,xposv,'b.',vdatafit,ydatafit,'black')
    axis([ymin*1.1 ymax*1.5 -framesize/2+(y_coord) framesize/2+(y_coord)])
    set(gca,'YDir','reverse')
    figure(gcf) % bring current figure to front
    hold off
    
    
end

end
function [] = HistPeakPlot(obj,frame_number)
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'

if nargin <2
    frame_number = 1;
end

for ii = frame_number;
    data = obj.Frame(:,:,ii);
    x_coord = obj.ApproxCentroid(1,1,ii);
    y_coord = obj.ApproxCentroid(1,2,ii);
    x(2) = obj.Centroid(1,1,ii); %xcen
    x(4) = obj.Centroid(1,2,ii); %ycen
    x(1) = obj.Centroid(1,3,ii); %peak
    x(3) = obj.Centroid(1,4,ii); %x width
    x(5) = obj.Centroid(1,5,ii); %y width
    x(6) = obj.Centroid(1,6,ii); %rotation
    
    
    framesize = round(obj.Centroid(1,4,ii)*10);
    
    [X,Y] = meshgrid(-framesize/2:framesize/2);
    X = X (1:framesize,1:framesize);
    Y = Y (1:framesize,1:framesize);
    
    X = X+x_coord;
    Y = Y+y_coord;
    
    Bot = Y(1,1);
    Top = Y(framesize,1);
    Left = X(1,1);
    Right = X(1,framesize);
    Z = frame(Bot:Top,Left:Right); % cut data at the locations corresponding to the chosen frame size.
    
    
    %figure creation
    hf2 = figure(ii+1);
    set(hf2, 'Position', [20 20 950 950])
    alpha(0)
    subplot(4,4,[5,6,7,9,10,11,13,14,15])
    imagesc(X(1,:),Y(:,1)',Z)
    set(gca,'YDir','reverse')
    colormap('bone')
    
    % -----Calculate cross sections-------------
    %generate points along horizontal axis
    m = -tan(x(6));% Point slope formula
    b = (-m*x(2) + x(4));
    xvh = -framesize/2+x_coord:framesize/2+x_coord;
    yvh = xvh*m + b;
    hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
    %generate points along vertical axis
    mrot = -m;
    brot = (mrot*x(4) - x(2));
    yvv = -framesize/2+y_coord:framesize/2+y_coord;
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
    axis([-framesize/2+(x_coord) framesize/2+(x_coord) -framesize/2+(y_coord) framesize/2+(y_coord)])
    
    %Histogram plots
    ymin = 0;
    
    xdatafit = linspace(-framesize/2+(x_coord), framesize/2+(x_coord),300);
    hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
    
    ydatafit = linspace(-framesize/2+(y_coord), framesize/2+(y_coord),300);
    vdatafit = x(1)*exp(-(ydatafit-x(4)).^2/(2*x(5)^2));
    
    subplot(4,4,[1:3])
    histogram(obj.Centroid(:,1),round(sqrt(length(obj.Centroid(:,1)))))
    xlim([obj.Mean(1,1)-1.5*(obj.Range(:,1)), obj.Mean(1,1)+1.5*(obj.Range(1,1))])
    
    subplot(4,4,[8,12,16])
    histogram(gca,obj.Centroid(:,2),round(sqrt(length(obj.Centroid(:,1)))))
    xlim([obj.Mean(1,2)-1.5*(obj.Range(1,2)), obj.Mean(1,2)+1.5*(obj.Range(1,2))])
    set(gca,'view',[90 -90])
    set(gca,'XDir','reverse')
    figure(gcf) % bring current figure to front
    hold off
end

end
function [] = ScatterPlot(obj)
figure()
hold on
axis equal
plot(obj.Centroid(:,1),obj.Centroid(:,2),'.')
box on
end
