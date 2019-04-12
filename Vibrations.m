%Vibrations

%--------------%
% SX, DX 
%--------------%
%(for frequency plot comparison)
load('S:\ImageProcessing\Results\New\SX2.mat');
load('S:\ImageProcessing\MLA\PSD_realtime.mat')
I = FGC;
c = 0.604*6.5;
% ----------------%
% Frequency
% ----------------%

% Fs = length(x)/2;
% 
% [f,psdx] = calcPSD(x,Fs);
% 
% [f,psdy] = calcPSD(y,Fs);

Fs = I.frameRate;

[f,psdx] = calcPSD(c*(I.fitParameters(:,2)-I.Mean(2)),Fs);
 
[f,psdy] = calcPSD(c*(I.fitParameters(:,4)-I.Mean(4)),Fs);

[f2,psdSR] = calcPSD(I.instSR',Fs);

csdy = 1e-3*sqrt(cumtrapz(f(2:end),psdy(2:end))); % arcseconds
csdy = csdy./max(csdy);

csdx = 1e-3*sqrt(cumtrapz(f(2:end),psdx(2:end))); % arcseconds
csdx = csdx./max(csdx);

figure()
subplot(1,2,1)
semilogy(f(2:end),(psdx(2:end)),'linewidth',1.5);
hold on
semilogy(f(2:end),(psdy(2:end)),'linewidth',1.5,'Color',[0.93,0.69,0.13]);
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
box on
grid on
xlabel('Frequency (Hz)')
ylabel('PSD mas^2/Hz')
xlim([min(f),200])
legend('Tip','Tilt');
legend boxoff

subplot(1,2,2)
plot(f(2:end),csdx,'linewidth',1.5);
hold on
plot(f(2:end),csdy,'linewidth',1.5,'Color',[0.93,0.69,0.13]);
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 1.5;
box on
grid on
xlabel('Frequency (Hz)')
ylabel('Cumulative PSD')
xlim([min(f),200])
legend('Tip','Tilt');
legend boxoff

function [f,psdx] = calcPSD(xdata,Fs)
N = length(xdata);
xdft = fft(xdata);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2; %alternative way (same answer) psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
f = Fs*(0:(N/2))/N; % alternative way (same answer) f2 = 0:Fs/N:Fs/2;
end