%efficiency time series 
clear
load('S:\ImageProcessing\Results\Andrew_reduced\NuCom4.mat')
load('S:\ImageProcessing\Results\Andrew_reduced\Aus153.mat')
load('S:\ImageProcessing\Results\Andrew_reduced\Aus30.mat')
load('S:\ImageProcessing\Results\Andrew_reduced\KXVir.mat')
load('S:\ImageProcessing\Results\Andrew_reduced\NSV.mat')


% Coupling efficiency plots 
science = viridis(6);

for ii = 1:length(science)-1
    colors{ii} = science(ii,:);
end

colors{ii+1} = [0.64, 0.08,0.18];
colors{ii+2} = [0, 0.45,0.74];
colors{ii+3} = [0.85, 0.33,0.1];

%Time series

t1 = (0:45000-1)*0.01; % 500 seconds
c1 = NuCom4.coupling(1:45000); 
e1(:,1) = (NuCom4.sigma).*ones(length(c1),1);

t2 = (0:38928-1)*0.01; % 
c2 = Aus153.coupling(1:38928); 
e2(:,1) = (Aus153.sigma).*ones(length(c2),1);

% t3 = (0:224998-1)*0.01; % 10Hz and number of data points
% c3 = Aus30.coupling(1:224998); 

t3 = (0:23503-1)*0.01; % reduce the length as this one drops from flexure quite a bit
c3 = Aus30.coupling(1:23503); 

e3(:,1) = (Aus30.sigma).*ones(length(c3),1);

t4 = (0:23503-1)*0.01; % 500 seconds
c4 = KXVir.coupling; 
e4(:,1) = (KXVir.sigma).*ones(length(c4),1);

t5 = (0:12743-1)*0.01; % 500 seconds
c5 = NSV.coupling; 
e5(:,1) = (NSV.sigma).*ones(length(c5),1);

% Enlarge figure to full screen.
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 1]);
hold on
% shadedErrorBar(t1,c1,e1,'lineprops',{'color',colors{1},'linewidth',1.5})
h(3) = shadedErrorBar(t2,100*c2,100*e2,'lineprops',{'color',colors{2},'linewidth',1.5});
h(2) = shadedErrorBar(t3,100*c3,100*e3,'lineprops',{'color',colors{3},'linewidth',1.5});
h(1) = shadedErrorBar(t4,100*c4,100*e4,'lineprops',{'color',colors{6},'linewidth',1.5});
% shadedErrorBar(t5,c5,e5,'lineprops',{'color',colors{7},'linewidth',1.5})
xlabel('Time (s)')
xlim ([0,120])
ylabel('Fiber Coupling (%)')
grid on
box on
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize = 22;
ylim([0,35])
l = legend([h(1).mainLine,h(2).mainLine,h(3).mainLine],{'400 modes HD 113496','300 modes HD 89758','153 modes HD 89758'});


%Frequency
% xdata = c2;
% Fs = 100;
% N = length(xdata);
% xdft = fft(xdata);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2; %alternative way (same answer) psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% f = Fs*(0:(N/2))/N; % alternative way (same answer) f2 = 0:Fs/N:Fs/2;
% 
% hf = figure();
% set(hf, 'Position', [20 20 1300 350])
% plot(f,psdx)
% xlabel('f (Hz)')
% ylabel('Power/ Frequency (Efficiency^2/Hz)')
% set(gca,'FontSize',16)
% xlim([0.01 1])



% % Enlarge figure to full screen.
% figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.625]);
% hold on
% shadedErrorBar(t3,c3,e3,'lineprops',{'color',[0.64,0.08,0.18],'linewidth',2})
% 
% xlabel('time(s)')
% % xlim ([0,360])
% ylabel('Fiber Coupling (%)')
% grid on
% box on
% ax = gca;
% ax.LineWidth = 1.5;
% ax.FontSize = 22;



