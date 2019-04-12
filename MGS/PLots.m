%MGS plots for PAPER

% Load saved figures
c=hgload('img2.fig');
k=hgload('opd_est.fig');

% Paste figures on the subplots

pause
A = allchild(get(c,'CurrentAxes'));
A = getimage(A);
figure()
h(1)=subplot(1,3,1);
imagesc(A);
axis image
box on
grid on
title('Measured PSF')
xlim([1187-50,1187+50])
ylim([933-50,933+50])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

pause
B = allchild(get(c,'CurrentAxes'));
B = getimage(B);
figure(3)
h(2)=subplot(1,3,2);
imagesc(B);
axis image
box on
grid on
title('Prediceted PSF')
xlim([512-50,512+50])
ylim([512-50,512+50])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

C = allchild(get(k,'CurrentAxes'));
C = getimage(C).*1e9;
figure(3)
h(3)=subplot(1,3,3);
imagesc(C)
axis image
box on
grid on
title('OPD rms = xx nm')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
c = colorbar('vert');
c.Label.String = 'WFE'
