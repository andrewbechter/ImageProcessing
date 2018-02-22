clear 
close

load('S:\AcquisitionCamera\SX_090717\focus1.mat')
f1 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus2.mat')
f2 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus3.mat')
f3 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus4.mat')
f4 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus5.mat')
f5 = I.Frame;
% load('S:\AcquisitionCamera\SX_090717\focus6.mat')
% f6 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus7.mat')
f7 = I.Frame;
% load('S:\AcquisitionCamera\SX_090717\focus8.mat')
% f8 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus9.mat')
f9 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus10.mat')
f10 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus11.mat')
f11 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus12.mat')
f12 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus13.mat')
f13 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus14.mat')
f14 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus15.mat')
f15 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus16.mat')
f16 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus17.mat')
f17 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus18.mat')
f18 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus19.mat')
f19 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus20.mat')
f20 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus21.mat')
f21 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus22.mat')
f22 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus23.mat')
f23 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus24.mat')
f24 = I.Frame;
load('S:\AcquisitionCamera\SX_090717\focus25.mat')
f25 = I.Frame;

total = f1+f2+f3+f4+f5+f7+f9+f10+f11+f12+f13+f14+f15+f16+f17+f18+f19+f20+f21+f22+f23+f24+f25;

figure
imagesc(total)
set(gca,'YDir','normal');
axis equal