% focus
clear 
close all

load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus1.mat')
f1 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus2.mat')
f2 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus3.mat')
f3 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus4.mat')
f4 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus5.mat')
f5 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus6.mat')
f6 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus7.mat')
f7 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus8.mat')
f8 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus9.mat')
f9 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus10.mat')
f10 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus11.mat')
f11 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus12.mat')
f12 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus13.mat')
f13 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus14.mat')
f14 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus15.mat')
f15 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus16.mat')
f16 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus17.mat')
f17 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus18.mat')
f18 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus19.mat')
f19 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus20.mat')
f20 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus21.mat')
f21 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus22.mat')
f22 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus23.mat')
f23 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus24.mat')
f24 = I;
load('S:\ImageProcessing\Reduced\Old_NIC\SX_090717\focus25.mat')
f25 = I;

% A = {f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25};
A = {f1,f2,f3,f4,f5,f6,f7,f8,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25};
figure
hold on
for ii = 1:20
    I = A{ii}; 
    x(ii) = I.Mean(1,1)*6.5e-3;
    y(ii) = I.Mean(1,2)*6.5e-3;
%     plot(x(ii),y(ii),'.')
end

x = x - x(1);
y = y - y(1);
z = [12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,11.5,10.5,9.5,8.5,7.5,6.5,5.5,4.5,3.5,2.5,1.5,0.5];
z(9:13)=[];

plot(x,z,'.')

total = f1.Frame./max(max(f1.Frame))+f2.Frame./max(max(f2.Frame))+f3.Frame./max(max(f3.Frame))...
    +f4.Frame./max(max(f4.Frame))+f5.Frame./max(max(f5.Frame))...
    +f7.Frame./max(max(f7.Frame))+f9.Frame./max(max(f9.Frame))...
    +f10.Frame./max(max(f10.Frame))+f11.Frame./max(max(f11.Frame))+f12.Frame./max(max(f12.Frame))...
    +f13.Frame./max(max(f13.Frame))+f14.Frame./max(max(f14.Frame))+f15.Frame./max(max(f15.Frame))...
    +f16.Frame./max(max(f16.Frame))+f17.Frame./max(max(f17.Frame))+f18.Frame./max(max(f18.Frame))...
    +f19.Frame./max(max(f19.Frame))+f20.Frame./max(max(f20.Frame))+f21.Frame./max(max(f21.Frame))...
    +f22.Frame./max(max(f22.Frame))+f23.Frame./max(max(f23.Frame))+f24.Frame./max(max(f24.Frame))...
    +f25.Frame./max(max(f25.Frame));
figure
imagesc(total)
set(gca,'YDir','normal');
axis equal