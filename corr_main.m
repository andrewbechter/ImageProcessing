%correlations
clear

%--------------%
% NSV (no overlap)
%--------------%
% load('S:\ImageProcessing\Results\New\NSV');
% NSV300 = NSV;
% load('S:\ImageProcessing\Results\Andrew_reduced\NSV');
% I = NSV300;
% F = NSV;

%--------------%
% Aus (no image data)
%--------------%
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus');
% F = Aus;

%--------------%
% Aus30 (1) (no time overlap)
%--------------%
% load('S:\ImageProcessing\Results\New\Aus1');
% I = Aus;
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus');
% F = Aus;

%--------------%
% Aus30 (1) (no overlap)
%--------------%
% load('S:\ImageProcessing\Results\New\Aus1');
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
% I = Aus1;
% F = Aus30;

% --------------%
% Aus30 (2)
% --------------%
% load('S:\ImageProcessing\Results\New\Aus2');
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
% I = Aus2;
% F = Aus30;

%--------------%
% Aust30 (3)
%--------------%
% load('S:\ImageProcessing\Results\New\Aus3');
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
% I = Aus3;
% F = Aus30;


%--------------%
% Aust30 (4)
%--------------%
% load('S:\ImageProcessing\Results\New\Aus4');
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
% I = Aus4;
% F = Aus30;


% --------------%
% Aus30 (5)
% --------------%
load('S:\ImageProcessing\Results\New\Aus5');
load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
I = Aus5;
F = Aus30;

%--------------%
% Aus30 (6) (no overlap)
%--------------%
% load('S:\ImageProcessing\Results\New\Aus6');
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus30');
% I = Aus6;
% F = Aus30;

%--------------%
% Aus153
%--------------%
% load('S:\ImageProcessing\Results\New\Aus153');
% Aus153Sub=Aus153;
% load('S:\ImageProcessing\Results\Andrew_reduced\Aus153');
% I = Aus153Sub;
% F = Aus153;

%--------------%
% NuCom300(no overlap)
%--------------%
% load('S:\ImageProcessing\Results\New\NuComSub300');
% load('S:\ImageProcessing\Results\Andrew_reduced\NuCom3');
% I = NuComSub300;
% F = NuCom3;

%--------------%
% NuCom400
%--------------%
% load('S:\ImageProcessing\Results\New\NuComSub400');
% load('S:\ImageProcessing\Results\Andrew_reduced\NuCom4');
% I = NuComSub400;
% F = NuCom4;

%--------------%
% KXVir
%--------------%

% load('S:\ImageProcessing\Results\New\KXVIRSub');
% load('S:\ImageProcessing\Results\Andrew_reduced\KXVir');
% KXVIRSub.time(3277:end,:) = [];
% KXVIRSub.instSR(:,3277:end) = [];
% KXVIRSub.fitParameters(3277:end,:) = [];
% I = KXVIRSub;
% F = KXVir;

%--------------------%
% PreProcess settings
%--------------------%

Fs = I.frameRate;

%Filter
% windowSize = floor(Fs/4); 
% s = (1/windowSize)*ones(1,windowSize);
% p = 1; % 1 seems to work well for normalized data...

%------------------------%
% Fiber data
%------------------------%

fiber_in(:,1) = F.time(:,1); % assign variables

fiber_in(:,2) = F.coupling;

%------------------------%
% Strehl ratio processing
%------------------------%
%smooths and builds time vectors (corrects fiber coupling sampling)

sr_in(:,1) = I.time(:,2); % assign variables

sr_in(:,2) = I.instSR';

% std(sr_in(:,2));

movs = Fs/(2*1);
% sr_in(:,2) = filter(s,p,sr_in(:,2)); % filters
sr_in(:,2) = movmean(sr_in(:,2),movs); % moving average

% std(sr_in(:,2));

[fiber_out,sr_out] = pre_process_time_series(fiber_in,sr_in,Fs); % (fiber data, ANDOR data, ANDOR sampling frq)

%---------------------%
% Tip-Tilt processing
%---------------------%
%smooths and built time vectors (corrects fiber coupling sampling)

%magnitude of combined tip-tilt
r = 2.33*sqrt((I.fitParameters(:,2)-I.Mean(2)).^2 + (I.fitParameters(:,4)-I.Mean(4)).^2);

tt_in(:,1) = I.time(:,2); % assign variables

tt_in(:,2) = r.*-1;

% tt_in(:,2) = filter(s,p,tt_in(:,2)); % filter
tt_in(:,2) = movmean(tt_in(:,2),movs); % moving average

[fiber1_out,tt_out] = pre_process_time_series(fiber_in,tt_in,Fs); % process


%just tip series

tt_in(:,2) = 2.33*(I.fitParameters(:,2)-I.Mean(2));%reassign

% tt_in(:,2) =filter(s,p,tt_in(:,2)); % filter  
tt_in(:,2) = movmean(tt_in(:,2),movs); % moving average

[~,tip] = pre_process_time_series(fiber_in,tt_in,Fs); %process 


%just tilt series
tt_in(:,2) = 2.33*(I.fitParameters(:,4)-I.Mean(4)); %reassign 

std(tt_in(:,2));

% tt_in(:,2) = filter(s,p,tt_in(:,2)); %filter 
tt_in(:,2) = movmean(tt_in(:,2),movs); % moving average

std(tt_in(:,2));

[~,tilt] = pre_process_time_series(fiber_in,tt_in,Fs); %process

% figure
% hold on
% plot(tip(:,2))
% plot(tilt(:,2))

%----------------%
% Cross Correlate
%----------------%

% Strehl ratio and Fiber (normalized)
[indDiff,acorr,lag,sig,f,a,t,stat] = correlate(fiber_out,sr_out,false,true);
[FC,SR,time,gamma] = final_conditioning(indDiff,acorr,lag,f,a,t,Fs,false);

% Strehl ratio and Fiber (non normalized)
[~,~,~,~,FC_series,SR_series,~,~] = correlate(fiber_out,sr_out,false,false);
[FC_series,SR_series,~,~] = final_conditioning(indDiff,acorr,lag,FC_series,SR_series,t,Fs,false);

% TT and fiber (normazlied)
[indDiff2,acorr2,lag2,sig2,f2,a2,t2,stat2] = correlate(fiber1_out,tt_out,false,true);
[FC2,TT,time2,gamma2] = final_conditioning(indDiff,acorr,lag,f2,a2,t2,Fs,false);

% Tip (non normalized)
[~,~,~,~,~,Tip_series,~,~] = correlate(fiber1_out,tip,false,false);
[~,Tip_series,~,~] = final_conditioning(indDiff,acorr,lag,f2,Tip_series,t2,Fs,false);

%Tilt (non normalized)
[~,~,~,~,~,Tilt_series,~,~] = correlate(fiber1_out,tilt,false,false);
[~,Tilt_series,~,~] = final_conditioning(indDiff,acorr,lag,f2,Tilt_series,t2,Fs,false);

%Combined TT and SR (filterd) with Fiber
comb = (TT+SR)./2;
gamma3 = corrcoef(FC,comb);
c = 2.33;

%--------------%
% SX, DX 
%--------------%
%(for frequency plot comparison)
% load('S:\ImageProcessing\Results\New\DX2.mat');
% I = DX2;
% c = 0.604*6.5;

% % ----------------%
% % Frequency
% % ----------------%
% 
% Fs = I.frameRate;
% 
% [f,psdx] = calcPSD(c*(I.fitParameters(:,2)-I.Mean(2)),Fs);
% 
% [f,psdy] = calcPSD(c*(I.fitParameters(:,4)-I.Mean(4)),Fs);
% 
% [f2,psdSR] = calcPSD(I.instSR',Fs);
% 
% csdy = 1e-3*sqrt(cumtrapz(f(2:end),psdy(2:end))); % arcseconds
% csdy = csdy./max(csdy);
% 
% csdx = 1e-3*sqrt(cumtrapz(f(2:end),psdx(2:end))); % arcseconds
% csdx = csdx./max(csdx);
% 
% figure()
% subplot(1,2,1)
% plot(f(2:end),psdx(2:end),'linewidth',1.5);
% hold on
% plot(f(2:end),psdy(2:end),'linewidth',1.5,'Color',[0.93,0.69,0.13]);
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.5;
% box on
% grid on
% xlabel('Frequency (Hz)')
% ylabel('PSD (mas^2/Hz)')
% xlim([min(f),65])
% legend('Tip','Tilt');
% legend boxoff
% 
% subplot(1,2,2)
% plot(f(2:end),csdx,'linewidth',1.5);
% hold on
% plot(f(2:end),csdy,'linewidth',1.5,'Color',[0.93,0.69,0.13]);
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.5;
% box on
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Cumulative PSD')
% xlim([min(f),65])
% legend('Tip','Tilt');
% legend boxoff


%-------------------%
% Plots Strehl (new)
%-------------------%

figure
plot(time,FC_series,'linewidth',1.2)
hold on
plot(time,SR_series,'linewidth',1.2)
xlabel('Time (s)')
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
xlim([min(time) max(time)])
ylim([0,0.30])


figure
subplot(2,1,1)
plot(time,SR_series,'linewidth',1.2)
xlabel('Time (s)')
ylabel('Strehl Ratio')
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
xlim([min(time) max(time)])
ylim([0,0.30])


subplot(2,1,2)
plot(time,SR,'linewidth',1.2)
xlabel('Time (s)')
ylabel('Normalized Units')
hold on
plot(time,FC,'linewidth',1.2)
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ylim(ax.YLim)
xlabel('Time (s)')
ylabel('Normalized Units')
legend('Strehl Ratio','Fiber Coupling')
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
xlim([min(time) max(time)])

%-------------------%
% Plots TT (new)
%-------------------%

figure
subplot(2,1,1)
hold on
plot(time,Tip_series,'linewidth',1)
plot(time,Tilt_series,'linewidth',1,'Color',[0.93,0.69,0.13])
xlabel('Time (s)')
ylabel('Tip-Tilt (mas)')
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
xlim([min(time) max(time)])
ylim([-20,20])
legend('Tip','Tilt')

subplot(2,1,2)
plot(time,TT,'linewidth',1.2)
xlabel('Time (s)')
ylabel('Normalized Units')
hold on
plot(time,FC,'linewidth',1.2)
box on
grid on
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ylim([-5,5])
xlabel('Time (s)')
ylabel('Normalized Units')
legend('Tip/Tilt','Fiber Coupling')
xlim([min(time) max(time)])


%----------------%
% Plots (original)
%----------------%
% 
% figure
% subplot(3,1,1)
% plot(time,FC)
% hold on
% plot(time,SR)
% xlabel('Time(s)')
% ylabel('Normalized Units')
% legend('Fiber Coupling','I-Band Strehl')
% box on
% grid on
% title(['Correlation =',num2str(100*gamma(1,2)),'%'])
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.5;
% 
% subplot(3,1,2)
% plot(time2,FC2)
% hold on
% plot(time2,TT)
% ylim(ax.YLim)
% xlabel('time(s)')
% ylabel('Normalized Units')
% legend('Fiber Coupling','Tip/Tilt')
% box on
% grid on
% title(['Correlation =',num2str(100*gamma2(1,2)),'%'])
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.5;
% 
% comb = (TT+SR)./2;
% gamma3 = corrcoef(FC,comb);
% 
% subplot(3,1,3)
% plot(time,FC2)
% hold on
% plot(time,comb)
% ylim(ax.YLim)
% xlabel('time(s)')
% ylabel('Normalized Units')
% legend('Fiber Coupling','SR+Tip/Tilt')
% box on
% grid on
% title(['Correlation =',num2str(100*gamma3(1,2)),'%'])
% ax = gca;
% ax.FontSize = 18;
% ax.LineWidth = 1.5;






function [f,psdx] = calcPSD(xdata,Fs)
N = length(xdata);
xdft = fft(xdata);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2; %alternative way (same answer) psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
f = Fs*(0:(N/2))/N; % alternative way (same answer) f2 = 0:Fs/N:Fs/2;
end