function [theta_x,theta_y]=TT_PSD_OIM(n,res)

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New description: 
% This function creates a time series of angular displacements for an OIM mirror to replicate the tip/tilt errors 
% expected from vibrations (or other effects) experienced at the LBT.
%
% INPUT: 
%  n  = array size
%  res = ??

% OUTPUT:
%  time - vector of time stamps 
%  theta_x - commanded motion of mirror in x-direction
%  theta_y - commanded motion of mirror in y-direction
%
% The program works by taking an input PSD and creating a randomized series of motions that well-mimics the deviations 
% expected in practice. We use the "clock time" provided by Matlab to send instructions dicipherable to the OIM mirror.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(mod(n,2)==1,'Error: array size must be odd.');
assert(mod(res,2)==1,'Error: res must be odd for symmetry.');
PSD=zeros(n);  % Initialize two-dimensional temporal PSD.
k0=3;          % Create knee in curve.
index=2;       % power law index
A=1;           % Arbitrary amplitude.
D=1;           % Angular deviation amplitude (DOUBLE CHECK UNITS OF THESE ...)
mid=(n-1)/2+1; % Mid-point of array.

rf=10;         % Resonant frequency (units of k0?).
B=6000*A;         % Resonant frequency amplitude (compare to A).
ind2=20;        % Index of power spectrum resonant frequency. Controls "sharpness" of resonant frequency. 
assert(mod(ind2,2)==0,'Error: invalid resonant frequency index.');

rf1=20;         % Resonant frequency (units of k0?).
B1=0;%5000*A;         % Resonant frequency amplitude (compare to A).
ind21=20;        % Index of power spectrum resonant frequency. Controls "sharpness" of resonant frequency. 
assert(mod(ind21,2)==0,'Error: invalid resonant frequency index.');

fprintf('\n   Generating power-spectral-density ... \n');
% Working in image plane
for h=1:n
    y=mid-h;
    for j=1:n        
        x=j-mid;
        k=2*sqrt(x^2+y^2)/res;           % Temporal frequency [cycles / time?].
        if (k*res/2 < n/4)               % Nyquist sample highest frequency. 
            PSD(h,j)=A/(1+(k/k0)^index)+B/(1+(k-rf)^ind2)+B1/(1+(k-rf1)^ind21); % Power-spectral-density with resonant frequency.
        end
    end
end

fprintf('   Converting to time series ... \n');
rand_phase=2*pi*rand(n);
surface=fft2(fftshift(sqrt(PSD).*exp(1i*rand_phase)));  % complex surface
surface0 = surface;
surface=imag(surface);  % surface error -- only need imaginary or real part!                                   

% Force mirror to behave with desired amplitude for angular deviations:
surface=D*surface/sqrt(sum(sum(surface.*surface))/(pi*(n/res/2)^2));

theta_x=surface(mid,:)'; % Define theta_x variations.
theta_y=surface(:,mid); % Define theta_y variations.

% Could save the file here. 

fprintf('Now checking that functions actually follow intended PSD ... \n');

figure(31)
% h=pcolor(sqrt(PSD));  
h = imagesc(sqrt(PSD));
colormap(gray(100));
%get(h);
%set(h,'LineStyle','none');
axis equal tight;
colorbar('vert');
title('Normalized Power Spectral Density','FontSize',14);
set(gca,'FontSize',14)

figure(41)
h=pcolor(fftshift(surface));  
colormap(gray(100));
get(h);
set(h,'LineStyle','none');
axis equal tight;
colorbar('vert');
title('Angular Deviation','FontSize',14);
set(gca,'FontSize',14)
xlabel('\theta_x [unitless]','FontSize',18);
ylabel('\theta_y [unitless]','FontSize',18);

figure(51)
plot(theta_x,'-b','LineWidth',2); % Define this as theta_x.
hold on;
plot(theta_y,'-k','LineWidth',2); % Define this as theta_y.
xlabel('Time [unitless]','FontSize',18);
ylabel('Angular Deviation [unitless]','FontSize',18);
legend('\theta_x','\theta_y');
set(gca,'FontSize',14)
xlim([1 n]);

%clear rand_phase PSD;

fprintf('   Lowest temporal frequency present: %-3.2f cycles per TIME (check units). \n',2/res);
fprintf('   Highest temporal frequency present: %-3.2f cycles per TIME (check units). \n',2*n/4/res);

% The following is true, and serves as a good reminder of how the code works. Key is that the # of cycles in the pupil plane 
% equals TWICE the pixel # in the PSD image plane. So, pixel #3 contributes a wave that repeats 6 times in the pupil plane 
% MATRIX (not aperture). The sampling can be derived from there.

fprintf('   Now checking whether results actually reproduce PSD ... \n');

% Pad array with zeros for FFT.
% tempx=zeros(1,n);
% tempy=zeros(1,n);
% indices=(mid-round(n/res/2):1:mid+round(n/res/2));
% tempx(indices)=theta_x(indices);
% tempy(indices)=theta_y(indices);

% Alternatively, can try to pad array directly instead of removing data.
tempx=padarray(theta_x,[0 (res-1)/2*n]);
tempy=padarray(theta_y,[0 (res-1)/2*n]);

% Take Fourier transform.
xpsd=fftshift(ifft(tempx));
ypsd=fftshift(ifft(tempy)); 

figure(61)
semilogy(PSD(mid,:),'-b','LineWidth',3)
hold on;
semilogy(sqrt(xpsd.*conj(xpsd)/sum(sum(xpsd.*conj(xpsd)))),'-r','LineWidth',1);
%plot(PSD(:,mid),'-k');
%plot([1:1:n],sqrt(ypsd.*conj(ypsd)/sum(sum(ypsd.*conj(ypsd)))),'--k','LineWidth',3);
axis tight
%legend('x_{PSD}','x_{gen}','y_{PSD}','y_{gen}');
set(gca,'FontSize',14)