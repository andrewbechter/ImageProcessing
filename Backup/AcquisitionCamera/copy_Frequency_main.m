clear
% close all

%% FFT Main (notes)
% This script handles reading in data files and taking the fourier
% transform. The plots are setup to show the time series data and single
% sided fourier transform. The limits are set to display half the sampling
% frequency (Nyquist sampling)

% Setting up the Time baselines and frequency sampling can be done either
% by using a known sampling frequency (Basler), known time intervals
% (Newport power meter) or from the actual time stamps and a
% recorded frame rate (ANDOR). The first two cases can be found in the main
% script, and the last case is more complicated but examples can be found
% in the ANDOR data sets.

% There is currently a low pass and high pass filter option (Low pass is
% moving average and high pass is an fir bandpass ripple. Only use the
% highpass filter if you understand the inputs required to make the correct
% filter for your data. The lowpass filter is simple to use. 

% A peak finding function and signal reconstruction section have yet to be
% added. 

%------------------------Data Sets------------------------%

%% Fiber Coupling (data format)
% Time, Time corrected, output power (W)
%% Tanaia Australis
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug18_Out_Australis_201015');
% data=Aug18_Out_Australis_201015(:,3);
% data= data(1:24000);
 
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug18_Out_Australis_212743');
% data=Aug18_Out_Australis_212743(:,3);
% data= data(1:150000);

%Small peaks (maybe due to small number of data points)
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug18_Out_Australis_220801');
% data=Aug18_Out_Australis_220801(:,3);
% data= data(1:12000);
%% Del_02
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug13_Out_Del02_040842');
% data=Aug13_Out_Del02_040842(:,3);
% data= data(1:12000); %loop opens after 15000, after 12000 there is weird doubling 

% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug13_Out_Del02_042131');
% data=Aug13_Out_Del02_042131(:,3);
% data= data(1:15000);%loop opens after 18000, after 15000 there is weird doubling 
%% Pollux
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Power_reduced\Aug14_Out_Pollux_201353');
% data=Aug14_Out_Pollux_201353(:,3);
% data= data(1:24000);

%% ANDOR (data format)
%Data Format:

%Column 1: Frame time microseconds
%Column 2: Actual time (based upon first frame timestamp to nearest second)
%Column 3: X pixel peak
%Column 4: Y pixel peak
%Column 5: Distance from median pixel based on Col 3 and Col 4
%Column 6: X sub-pixel peak (Gaussian fit)
%Column 7: Y sub-pixel peak (Gaussian fit)
%Column 8: Distance from median pixel based on Col 6 and Col 7 
%Column 9: theta, rotation angle in radians 
%Column 10: max_count  peak count of the Gaussian fit 
%Column 11: x centoid guess (variable used in script is x_coord) 
%Column 12: y centroid guess (variable used in script is y_coord)
%Column 13: x_sigma, sigma is the width parameter (2.35*sigma = FHWM) 
%Column 14: y_sigma 
%Column 15: Bad frame
%% Tania Australis
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\ANDOR\ANDOR_reduced\Apr18_Australis_2')
% data = Apr18_Australis_2(:,7);
% t = Apr18_Australis_2(:,1);
% t(:,2) = Apr18_Australis_2(:,2);
% L = length(data);      % Length of signal
% time(:,1) = t(:,1)*1e-6;
% Fs = 123.78;
%% KxVir
%load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\ANDOR\ANDOR_reduced\Aug18_KXVir_Sub')
% data = Apr18_KXVir_Sub(1:2000,6);
% t = Apr18_KXVir_Sub(1:2000,1);
% t(:,2) = Apr18_KXVir_Sub(1:2000,2);
% L = length(data);      % Length of signal
% time(:,1) = t(:,1)*1e-6;
% Fs = 58.47;

%% NuCom
% load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\ANDOR\ANDOR_reduced\Apr18_NuCom_Sub_400')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\ANDOR\ANDOR_reduced\6-20-17_Andrew\Apr18_NuCom_Sub_400')
data = Apr18_NuCom_Sub_400(:,6);
% t = Apr18_NuCom_Sub_400(:,1);
% t(:,2) = Apr18_NuCom_Sub_400(:,2);
% L = length(data);      % Length of signal
% time(:,1) = t(:,1)*1e-6;
Fs = 130;  % Sampling frequency
T = 1/Fs;             % Sampling period
L = length(data);      % Length of signal
time(:,1) = (0:L-1)*T;

%% Lab (test data)
%% Noise
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\Noise');
% data=power(:,1);
%% No Ocillation 
% high power
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\No_oscillation_high');
% data=power(:,1);

% star power
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\No_oscillation_low');
% data=power(:,1);
%% 10Hz Ocillation high power high amplitude (new code 12ms sampling???)
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\10Hz_high_007');
% data=power(:,1);

%10Hz Ocillation high power low amplitude 
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\10Hz_high_001');
% data=power(:,1);

%10Hz Ocillation high power highest amplitude 
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\10Hz_high_014');
% data=power(:,1);
%% 5Hz Ocillation (original code 10ms sampling)
% high power high amplitude 
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\5Hz_high_007');
% data=power(:,1);
%% 10Hz Ocillation (original code 10ms sampling)
%high power high amplitude
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\10Hz_high_007_orig');
% data=power(:,1);

%low power high amplitude
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\10Hz_low_007_orig');
% data=power(:,1);

%star power high amplitude
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\PowerMeter\10Hz_star_007_orig');
% data=power(:,1);
%% Basler
%%10Hz Ocsillation (picoscope basler 90um)
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\Basler\basler_10hz_50hz');
% data=basler_10hz_50hz(:,:,1);

% Timed basler recording 20ms = 50hz. 
% load('S:\Demonstrator_OnSky_Analysis\Frequency\Lab_FFT\Basler\10hz_50fps_timed');
% data=h(:,:,1);

% Timed basler recording
% load('P:\iLocater\iLocater_Demonstrator\Lab_Stability\Set16\basler_120ms_30s');
% data=basler_120ms_30s(:,1);

% Timed basler recording
% load('P:\iLocater\iLocater_Demonstrator\Lab_Stability\Set18\basler_10ms_10s');
% datah=basler_10ms_10s(:,1);
% mean_h = mean(datah);
% disph=datah-mean_h;
% datav=basler_10ms_10s(:,2);
% mean_v = mean(datav);
% dispv=datav-mean_v;
% data = sqrt(abs(disph.^2)+abs(dispv.^2));

%% ------------ Time Sampling Information ------------------------- %%
%% Known sampling frequency
% Fs = 50;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = length(data);      % Length of signal
% time(:,1) = (0:L-1)*T;

%% Known time intervals
% T = 0.01;
% Fs = 1/T;
% L = length(data);      % Length of signal
% time(:,1) = (0:L-1)*T;
%% ----------------------- AO Data ----------------------------%

%AO_name = 'AO_18_slope';
% AO_name = 'AO_18_ttres';
% AO_file = strcat('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\INDI\AO_reduced\',AO_name,'.mat');
% load(AO_file,AO_name);
% AO_data = eval(AO_name);
% 
% if(AO_data(1,1) <= t(1,2))
%     AO_ind_s = 1;
%     AO_time = AO_data(AO_ind_s,1);
%         while(AO_time <= t(1,2))
%             AO_ind_s = AO_ind_s + 1;
%             AO_time = AO_data(AO_ind_s,1);
%             if(mod(AO_ind_s,100000) == 0)
%                          fprintf('%i\n',AO_ind_s)
%             end
%         end
% else
%     AO_start_s = 1;
% end
% 
% if(AO_data(end,1) >= t(end,2))
%     AO_ind_e = size(AO_data,1);
%     AO_time = AO_data(AO_ind_e,1);
%         while(AO_time <= t(1,2))
%             AO_ind_e = AO_ind_e - 1;
%             AO_time = AO_data(AO_ind_e,1);
%             if(mod(AO_ind_e,100000) == 0)
%                          fprintf('%i\n',AO_ind_e)
%             end
%         end
% else
%     AO_start_e = size(AO_data,1);
% end
% 
% %data = AO_data(AO_ind_s:AO_ind_e,2);
% %t = (AO_data(AO_ind_s:AO_ind_e,1) - AO_data(AO_ind_s,1)) * (24 * 60 * 60);
% 
% data = AO_data(AO_ind_s+1000:AO_ind_s+40000,2);
% t = (AO_data(AO_ind_s+1000:AO_ind_s+40000,1) - AO_data(AO_ind_s,1)) * (24 * 60 * 60);
% 
% L = length(data);      % Length of signal
% time = t(:,1);
% T = time(2,1)-time(1,1);
% Fs = 1/T;

%% ------------------------Main Script------------------------%
%% Filters
% Moving Average
%windowSize = 6;
%[smooth_data] = Moving_Average(windowSize,data);

% High Pass Fir Filter
% F_Stop is stop band frequency
% Fpass is passband frequency
% Astop is stopband Attenuation
% Apass is PassbandRipple

% p = [Fstop, Fpass, Astop, Apass, Sampling] 
% p = [1,2,1,1,Fs];
% [high_data] = high_pass_fir(p,data);

%% Fourier transform of filtered or unfiltered data
[f,P1]=fourier_transform (data,L,Fs); %output = [frequency, powerspectrum]

return
%% Plots
figure()
plot(time(:,1),data)

figure()
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0.5 Fs/2])


