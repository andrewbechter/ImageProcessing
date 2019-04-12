% Coupling efficiency plots 
science = viridis(8);

for ii = 1:length(science)-1
    colors{ii} = science(ii,:);
end


load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\Aus.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\Aus153.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\KXVir.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\Leo72.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\NSV.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\NuCom3.mat')
load('P:\iLocater\iLocater_Demonstrator\LBT_Data\Matlab\Power_Data\Andrew_reduced\NuCom4.mat')


%x = [14.3,17.6,12,24.1,28.9,27.9,7];% old vales
norm = 0.78;
x = [19.2,19.5,25.7,26.8,26,12,7];% new values
y = [Aus153.eff,Aus.eff,NuCom3.eff,NuCom4.eff,KXVir.eff, Leo72.eff,NSV.eff]*100;
e = [Aus153.sigma,Aus.sigma,NuCom3.sigma,NuCom4.sigma,KXVir.sigma, Leo72.sigma,NSV.sigma]*100;

lambda = 800e-9;
phi = sqrt(abs(log(x/100))).*lambda./1020e-9;
strehl_ratio = 100*exp(-(phi.^2));

% figure
% ax1=gca ;
% errorbar(x,y,e,'.','markersize',15,'linewidth',1.5,'Parent',ax1,'color',[0.64,0.08,0.18])
% hold on
% grid on
% xlabel('I-band Strehl Ratio (%)')
% ylabel('Y-band SMF coupling efficiency (%)')

figure
hold on
ax1=gca;
for ii = 1:length(strehl_ratio)
h(ii)=errorbar(strehl_ratio(ii),y(ii)./norm,e(ii),'.','markersize',15,'linewidth',1.5,'Parent',ax1,'color',[0.64,0.08,0.18]);
end
% xlim([15,45])
% ylim([15,45])
hold on
grid on
box on
ax1.FontSize = 16;
xlim([15,50])
ft = fittype('a*x');
f = fit(strehl_ratio',(y./norm)',ft);
l2 = plot(f);
l2.Color = [0.64,0.08,0.18];
l2.LineWidth = [1.5];
l2.LineStyle = '--';
xlabel('Y-band Strehl Ratio (%)')
ylabel('Norm. coupling efficiency (%)')
ax1.LineWidth = 1;
label(l2,['slope = ',num2str(f(1))],'slope','center')
c = {'HD 89758 (153)','HD 89758 (300)','HD 113496 (300)','HD 113496 (400)','SAO 82686','HD 97778','NSV 19434'};
dx = 1; dy = 0; % displacement so the text does not overlay the data points
text(strehl_ratio+dx, y./norm+dy, c,'FontSize',14');
% legend([h(1),h(2),h(3),h(4),h(5),h(6),h(7)],c)

% 4 axis setup

% figure
% ax1=gca ;
% ax1_pos = ax1.Position; % position of first axes
% errorbar(x,y,e,'.','markersize',15,'linewidth',1.5,'Parent',ax1,'color',[0.64,0.08,0.18])
% hold on
% grid on
% xlabel('I-band Strehl Ratio (%)')
% ylabel('Y-band SMF coupling efficiency (%)')
% 
% % setup second x axis
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% lb = min(ax1.XLim);
% lb = sqrt(abs(log(lb/100))).*lambda./1020e-9;
% lb = 100*exp(-(lb.^2));
% ub = max(ax1.XLim);
% ub = sqrt(abs(log(ub/100))).*lambda./1020e-9;
% ub = 100*exp(-(ub.^2));
% ax2.XLim = ([lb,ub]);
% ax2.XLabel.String = ('Y-band Strehl Ratio (%)');
% % setup second y axis
% 
% lb = min(ax1.YLim);
% lb = lb./norm;
% ub = max(ax1.YLim);
% ub = ub./norm;
% ax2.YLim = ([lb,ub]);
% ax2.YLabel.String = ('Normalized Y-band coupling (%)');
% ax = ax2;