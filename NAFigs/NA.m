clear all
close
% Numerical aperture: Calculates the numerical aperture from measured data

%% Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end
colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
clear d

figure(1)

%% Old Data path shit

% path950_1 = 'P:\iLocater\GSFC_Fibers\NA\950nm\8_31\Sig1\';
% path950_2 = 'P:\iLocater\GSFC_Fibers\NA\950nm\8_31\Sig2\';
% 
% path980_1 = 'P:\iLocater\GSFC_Fibers\NA\980nm\8_31\Sig1\';
% path980_2 = 'P:\iLocater\GSFC_Fibers\NA\980nm\8_31\Sig2\';
% path980_3 = 'P:\iLocater\GSFC_Fibers\NA\980nm\8_31\Sig3\';
% 
% path1050_1 = 'P:\iLocater\GSFC_Fibers\NA\1050nm\8_31\Sig1\';
% path1050_2 = 'P:\iLocater\GSFC_Fibers\NA\1050nm\8_31\Sig2\';
% path1050_3 = 'P:\iLocater\GSFC_Fibers\NA\1050nm\8_31\Sig3\';
% 
% path950= 'P:\iLocater\GSFC_Fibers\NA\950nm\8_29\';
% path980 = 'P:\iLocater\GSFC_Fibers\NA\980nm\8_29\';
% path1050 = 'P:\iLocater\GSFC_Fibers\NA\1050nm\8_29\';
% path1064 = 'P:\iLocater\GSFC_Fibers\NA\1064nm\8_29\';
%%
path950_new = 'P:\iLocater\GSFC_Fibers\NA2\950nm\9_1\';
path980_new = 'P:\iLocater\GSFC_Fibers\NA2\980nm\9_1\';
path1050_new = 'P:\iLocater\GSFC_Fibers\NA2\1050nm\9_1\';

ii = 1;
% for path_list = {path950_new,path980_new,path1050_new,path950,path980,...
% path1050,path950_1,path980_1,path1050_1}%path1064}

for path_list = {path950_new,path980_new,path1050_new}

path = path_list{1};   
load([path,'0mm']);load([path,'10mm']);load([path,'20mm']);load([path,'30mm']);load([path,'40mm']);
pix = 5.4;
sigma = [];
width = [];
error = [];
RHSwidth =[];

for list = {I_0,I_10,I_20,I_30,I_40}
    dummy = abs(list{1}.Centroid(1,5)*pix); %Gaussian method
    sigma = [sigma dummy]; %Gaussian method
    
    e2 = (list{1}.Frame(round(list{1}.Centroid(2)),:)'); % use a slice
    e2 = mean(list{1}.Frame,1)'; % use the mean profile
    e2 = mean(list{1}.Frame,2); % use the mean profile
    e2=e2./(max(e2));%normalize the profile
    x = 1:length(e2);%pixel values
    ref = 1/exp(2);%reference point we are looking for on the profile
    ref = 0.06;
    % Fitting and finding the ref intesnity% values
    [a,b,c]=PVoigt((1:length(e2)),e2);   
    x2 = x-a(2);%centered x values
    
    
    %LHS width (cant always use LHS because wide images fall off LHS edge.
%     [ind,val]= find(c(x)>(ref),1,'first');
%     ind2 = ind-1;
%     LHSref = interp1([c(ind2),c(ind)],[ind2,ind],ref);
%     dummy2 = abs(a(2)-LHSref)*pix;
%     width = [width dummy2]; 
%     rms = sqrt(mean((e2-c(x)).^2)/2);% fit error 
       
    %LHS width error  
%     [ind4,val3] = find(c(x)>(ref+rms),1,'last'); 
%     ind5 = ind4+1;
%     rhs_error = interp1([c(ind4),c(ind5)],[ind4,ind5],ref+rms);
%     width_error =abs(rhs_error-a(2));
%     fit_error = abs(abs(a(2)-LHSref)-width_error);
    

    x3 = (1:4000);%extended x range while keeping the center point valid
    %RHS width
    [Rind,val]= find(c(x3)>(ref),1,'last');
    Rind2 = Rind+1;
    RHSref = interp1([c(Rind),c(Rind2)],[Rind,Rind2],ref);
    dummy3 = abs(a(2)-RHSref)*pix;
    width = [width dummy3]; 
    rms = sqrt(mean((e2-c(x)).^2)/2);% fit error 
       
    %RHS width error (cant use LHS because wide images fall off LHS edge. 
    [Rind3,val3] = find(c(x3)>(ref+rms),1,'last'); 
    Rind4 = Rind3+1;
    rhs_error = interp1([c(Rind3),c(Rind4)],[Rind3,Rind4],ref+rms);
    width_error =abs(rhs_error-a(2));   
    fit_error = abs(abs(a(2)-RHSref)-width_error);
    
    %Errors per data point
    error_exp = sqrt(0.1^2 + 0.5^2); %pixel size error and postioning error (tranlated to w in microns)
    fit_error = fit_error*pix; % fit error in microns
    ph_error = []; %how can I asses the photon noise error?
    dummy_error = sqrt(error_exp^2+fit_error^2); %total error in each data point
    error = [error dummy_error];
    
    figure(2)
    hold on
    plot(x2,e2,'.')
    plot(x2,c(x))
    xlabel('Distance (pixels)')
    ylabel('Normalized Intensity')
    box on
    grid on
    %% Bessel function stuff   
%     x3 = 2*(x-a(2))./length(x);
%     keyboard
%     e3 = e2(ind:ind3);
%     x4 = linspace(-1,1,length(e3));
%     %create fitting function and options
%     fitFunc = @(p,x) besselj(0,p*x).^2;
%     options = fitoptions(...
%     'Method', 'NonlinearLeastSquares',...
%     'TolX',1e-7);
% 
%     %fit the data and get parameter
%     [fitobject,gof,output] = fit(x4',e3,fitFunc,options)
%     parameter = fitobject.p;
%     confidence_interval = confint(fitobject);
% 
%     %plot
%     figure
%     plot(x4',e3,'.');
%     hold on;
%     plot(x4,besselj(0,1.2*x4).^2)
%     plot(x4,fitobject(x4),'k-','LineWidth',1.3);
%     set(gca,'FontSize',14);
%     title(sprintf('Fit with p = %.2f  [%.2f, %.2f]',parameter, confidence_interval));
%     legend('data','fit');     
    
end
figure(1)
hold on
z = [0,10,20,30,40]*1e3; % distance between sample points in microns
w = width-width(1);% take the difference in the width from the starting point
lin= 'a1*x';%linear model with origin as y intercept
weight = (1./error).^2;
[f,gof]=fit(z',w',lin,'Start',[0],'Weights',weight); % fitting the function
y = f(z');%model y-data
errorbarxy(z,w,[],error,[],error,'.',colors{ii}) % plot error bars 

hold on
h(ii)=plot(z,y,'Color',colors{ii}); %plots with figure handles for legened
p = coeffvalues(f); % slope coefficieny
e = confint(f);%condifence interval of the NA
na(ii,1) = p;
na(ii,2) = mean(abs(e-p));
ii = ii+1;
end

xlabel('Distance (\mum)')
ylabel('PSF radius (2\sigma) in \mum')
box on
grid on
% l=legend(h,{'950nm','980nm','1050nm','950nm','980nm','1050nm','950nm 1','980nm 1','1050nm 1'});%,'1064nm'});
l=legend(h,{'950nm','980nm','1050nm'});%,'1064nm'});
ax = gca;
ax.GridAlpha = 0.1;
ax.GridLineStyle ='--';
ax.LineWidth = 0.75;
l.Color = 'none';
l.Location = 'Best';

% Create textbox
annotation(figure(1),'textbox',[0.22 0.7 0.35 0.145],...
    'String',{['NA =',num2str(na(1,1)),' @\lambda = 950nm'],['NA =', num2str(na(2,1)),' @\lambda = 980nm'],['NA =', num2str(na(3,1)),' @\lambda = 1050nm']});



