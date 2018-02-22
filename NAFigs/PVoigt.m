function [params,good,f] = PVoigt(xdata,ydata)   


startpoints = [0 0 0 0 0];

    startpoints(1) = max(ydata);
    startpoints(2) = length(xdata)/2;
    startpoints(3) = 200; 
    startpoints(4) = 0.5;

    fitfunc = 'a * (d * exp(-0.5 * (x-b).^2 ./ (c/(2*sqrt(2*log(2)))).^2)+(1-d)*(c/2)^2./((x-b).^2+ (c/2)^2))+e';%*x.^2 +f*x +g';
    
    
    fo = fitoptions('Method','NonlinearLeastSquares','Start',startpoints, 'Lower', [0 0 0 0 0],'Upper', [inf inf inf 1 inf]);%,'display','notify','TolX',1e-6,'TolFun',1e-6);%,'MaxIter',600);%,'display','final');
    ft = fittype(fitfunc,'options',fo);
    [myfit,gof] = fit(xdata',ydata,ft);

    good = gof.adjrsquare;
    params = coeffvalues(myfit);
    f = myfit;
end 