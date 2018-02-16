clear all
% Numerical aperture: Calculates the numerical aperture from measured data

% path = 'P:\iLocater\GSFC_Fibers\NA\HighRes\';
% load([path,'1mm']); load([path,'3mm']);load([path,'5mm']);load([path,'7mm']);
% load([path,'9mm']);load([path,'11mm']);load([path,'13mm']);load([path,'15mm']);
% load([path,'17mm']);load([path,'19mm']);load([path,'21mm']);load([path,'23mm']);
% load([path,'25mm']);
% 
% pix = 5.4e-6;
% sigma = [];
% 
% for list = {I_1,I_3,I_5,I_7,I_9,I_11,I_13,I_15,I_17,I_19,I_21,I_23,I_25};
%     dummy = abs(list{1}.Centroid(1,4)*pix);
%     sigma = [sigma dummy];
% end
% 
% clear I_1 I_3 I_5 I_7 I_9 I_11 I_13 I_15 I_17 I_19 I_21 I_23 I_25
% 
% d = (1:2:25)*1e-3;
% delta_sigma = diff(sigma);
% Dist = d-d(1);
% Height =2*(sigma-sigma(1));
% 
% figure
% plot(Dist,Height,'.')

% path = 'P:\iLocater\GSFC_Fibers\NA\HighRes2\';
% load([path,'0mm']); load([path,'5mm']);load([path,'10mm']);load([path,'15mm']);load([path,'20mm']);
% load([path,'25mm']);
% 
% pix = 5.4e-6;
% sigma = [];
% for list = {I_0,I_5,I_10,I_15,I_20,I_25};
%     dummy = abs(list{1}.Centroid(1,4)*pix);
%     sigma = [sigma dummy];
% end
% 
% offset = (13+14)*1e-3;
% Dist = (0:1:5)*5e-3+offset;
% 
% delta_sigma = diff(sigma);
% % Height =2*(sigma-sigma(1));
% Height =2*sigma;
% 
% figure
% plot(Dist,Height,'.')


path = 'P:\iLocater\GSFC_Fibers\NA\HighRes3\';
load([path,'0mm']);load([path,'10mm']);load([path,'20mm']);

pix = 5.4e-6;
sigma = [];
for list = {I_0,I_10,I_20};
    dummy = abs(list{1}.Centroid(1,4)*pix);
    sigma = [sigma dummy];
end

Dist = [0,10,20]*1e-3;
delta_sigma = diff(sigma);
Height =2*(sigma-sigma(1));

figure
plot(Dist,Height,'.')