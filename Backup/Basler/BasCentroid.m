function [x_coord,y_coord,x_centroid,y_centroid,sigma_x,y_sigma,max_count,theta] = BasCentroid (frame,M,jj)

FitForOrientation=1;
% imagesc(data)
% data=double(data_unit8);
data=frame;
data_real = data;

%%%%%%%%%%%%%%%%%%%%---------------Start of gaussian fitting-----------------%%%%%%%%%%%%%%%%%%%%%%%
peakamp = max(max(data));
[C,I] = max(data(:));
C;
data(I);
[I1,I2] = ind2sub(size(data),I);
data(I1,I2);
x_coord = I2;
y_coord = I1;

samp = 0.1;
bfindbr= [floor(size(data,2)*(1-samp)), size(data,2),floor(size(data,1)*(1-samp)), size(data,1)]; 
bfindtr= [floor(size(data,2)*(1-samp)), size(data,2),1,floor(samp*size(data,1))]; 
bfindbl= [1, floor(samp*size(data,2)),floor(size(data,1)*(1-samp)), size(data,1)] ;
bfindtl= [1, floor(samp*size(data,2)),1,floor(samp*size(data,1))] ;


br = mean(mean(data(bfindbr(3):bfindbr(4),bfindbr(1):bfindbr(2))));
tr = mean(mean(data(bfindtr(3):bfindtr(4),bfindtr(1):bfindtr(2))));
bl = mean(mean(data(bfindbl(3):bfindbl(4),bfindbl(1):bfindbl(2))));
tl = mean(mean(data(bfindtl(3):bfindtl(4),bfindtl(1):bfindtl(2))));


if y_coord > (size(data,1)/2)
    if x_coord>(size(data,2)/2)
        %disp('botright')
       noise= mean([tr,bl,tl]);
    else
        %disp('botleft')
       noise=mean([tr,br,tl]);
    end
    
else 
    if x_coord>(size(data,2)/2);
        %disp('topleft')
       noise= mean([tr,br,bl]);
    else
        %disp('topright')
       noise= mean([tl,br,bl]);
    end
end

%%%%%%%%%%%----------------------noise reduction----------------------%%%%%%%%%%%
data = data-noise;

framesize=400;
[X,Y] = meshgrid(-framesize/2:framesize/2);
X = X (1:framesize,1:framesize);
Y = Y (1:framesize,1:framesize);

X = X+x_coord;
Y = Y+y_coord;

xdata = zeros(size(X,1),size(Y,2),2);

xdata(:,:,1) = X;
xdata(:,:,2) = Y;

Bot = x_coord-(framesize/2)+1;
Top = x_coord+(framesize/2);
Left = y_coord-(framesize/2)+1;
Right = y_coord+(framesize/2);

Z = data(Left:Right,Bot:Top);
Z=double(Z);

x0 = [peakamp,x_coord,9,y_coord,9,0]; %Inital guess parameters stored in array x0
xin = x0; 

InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
options=optimset('Diagnostics','off','Display','none');

% --- Fit---------------------
if FitForOrientation == 1
    % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
    lb = [peakamp,x_coord-5,1,y_coord-5,1,-pi/4];
    ub = [peakamp*2,x_coord+5,100,y_coord+5,100,+pi/4];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@BasD2GaussFunctionRot,x0,xdata,Z,lb,ub,options);
else
    x0 =x0(1:5);
    x0(6) = 0; 
    lb = [0,x_coord-5,1,y_coord-5,1,-inf];
    ub = [peakamp*10,x_coord+5,1000,y_coord+5,1000,+inf];
    [x,resnorm,residual,exitflag] = lsqcurvefit(@BasD2GaussFunction,x0,xdata,Z,lb,ub,options);
    %x(6) = 0;
end

dispangle = x(6)*180/pi;


% x(1)=193;
% x(2)=633.921;
% x(3)=1.36;
% x(4)=1256.424;
% x(5)=1.36;

% x0 =x0(1:6);
% lb = [peakamp*0.9,x_coord-2,1,y_coord-2,1,-pi/4];
% ub = [peakamp*1.1,x_coord+2,20,y_coord+2,15,+pi/4];
% %[x,resnorm,residual] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub,options);
% [x,resnorm,residual] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub,options); 



%%%%%%%%%%--------------JUST FOR PLOTTING-------------------%%%%%%%%%%%%%%%%%%%

if jj ==M

framesize=100;
[X,Y] = meshgrid(-framesize/2:framesize/2);
X = X (1:framesize,1:framesize);
Y = Y (1:framesize,1:framesize);

X = X+x_coord;
Y = Y+y_coord;

Bot = x_coord-framesize/2+1;
Top = x_coord+framesize/2;
Left = y_coord-framesize/2+1;
Right = y_coord+framesize/2;
Z = data(Left:Right,Bot:Top);

MdataSize=framesize;

hf2 = figure(1);
set(hf2, 'Position', [20 20 950 900])
alpha(0)
subplot(4,4, [5,6,7,9,10,11,13,14,15])
imagesc(X(1,:),Y(:,1)',Z)
set(gca,'YDir','reverse')
colormap('jet')

string1 = ['       Amplitude','    X-Coordinate', '    X-Width','    Y-Coordinate','    Y-Width','     Angle'];
string2 = ['Set  ',num2str(xin(1),'% 100.3f'),'           ',num2str(xin(2),'% 100.3f'),'       ',num2str(xin(3),'% 100.3f'),'         ',num2str(xin(4), '% 100.3f'),'        ',num2str(xin(5), '% 100.3f'),'     ',num2str(dispangle, '% 100.3f')];
string3 = ['Fit   ',num2str(x(1),'% 100.3f'),'            ',num2str(x(2),'% 100.3f'),'      ',num2str(x(3),'% 100.3f'),'         ',num2str(x(4), '% 100.3f'),'        ',num2str(x(5), '% 100.3f'),'     ',num2str(dispangle, '% 100.3f')];

text((x_coord)-MdataSize/3,+((y_coord)+MdataSize/2)+25,string1,'Color','k')
text((x_coord)-MdataSize/3,+((y_coord)+MdataSize/2)+45,string2,'Color','red')
text((x_coord)-MdataSize/3,+((y_coord)+MdataSize/2)+65,string3,'Color','red')

% -----Calculate cross sections-------------
%generate points along horizontal axis
m = -tan(x(6));% Point slope formula
b = (-m*x(2) + x(4));
xvh = -MdataSize/2+x_coord:MdataSize/2+x_coord;
yvh = xvh*m + b;
hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
%generate points along vertical axis
mrot = -m;
brot = (mrot*x(4) - x(2));
yvv = -MdataSize/2+y_coord:MdataSize/2+y_coord;
xvv = yvv*mrot - brot;
vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);

hold on % Indicate major and minor axis on plot
% 
% % plot points 
plot(xvh,yvh,'r.') 
plot(xvv,yvv,'g.')
% 
%plot lines 
plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'color',[0 .4 0]) 
% 
hold off
axis([-MdataSize/2+(x_coord) MdataSize/2+(x_coord) -MdataSize/2+(y_coord) MdataSize/2+(y_coord)])
%

ymin = 0;
ymax = x(1);

xdatafit = linspace(-MdataSize/2+(x_coord), MdataSize/2+(x_coord),300);
hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));

ydatafit = linspace(-MdataSize/2+(y_coord), MdataSize/2+(y_coord),300);
vdatafit = x(1)*exp(-(ydatafit-x(4)).^2/(2*x(5)^2));

subplot(4,4,[1:3])
xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0

plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
axis([-MdataSize/2+(x_coord) MdataSize/2+(x_coord) ymin*1.1 ymax*1.5])
subplot(4,4,[8,12,16])
xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0

plot(vPoints,xposv,'b.',vdatafit,ydatafit,'black')
axis([ymin*1.1 ymax*1.5 -MdataSize/2+(y_coord) MdataSize/2+(y_coord)])
set(gca,'YDir','reverse')
figure(gcf) % bring current figure to front
hold off


else
end


max_count= x(1);
x_centroid= x(2);
sigma_x = x(3);
y_centroid= x(4);
y_sigma= x(5);
theta= x(6);


end






