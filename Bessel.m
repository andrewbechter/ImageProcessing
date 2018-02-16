%for AngryBlackMan33
%in this example we fit the function besseli(1,p*x) where p=3
%create data that we want to fit
xs = linspace(-10,10)'; %apostrophe necessary: xs and data have to be vectors with shape (n,1)
data = besselj(0,1*xs); %data is a besseli function and some random noise
%plot
figure('Color','white');
plot(xs,data,'bo');
return
%create fitting function and options
fitFunc = @(p,x) besselj(0,p*x);
options = fitoptions(...
    'Method', 'NonlinearLeastSquares',...
    'Start',5,...
    'Lower',0,...
    'Upper',10,...
    'TolX',1e-7);

%fit the data and get parameter
[fitobject,gof,output] = fit(xs,data,fitFunc,options)
parameter = fitobject.p;
confidence_interval = confint(fitobject);

%plot
figure('Color','white');
plot(xs,data,'bo');
hold on;
plot(xs,fitobject(xs),'k-','LineWidth',1.3);
set(gca,'FontSize',14);
title(sprintf('Fit with p = %.2f  [%.2f, %.2f]',parameter, confidence_interval));
legend('data','fit');