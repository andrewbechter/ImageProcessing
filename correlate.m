function [lagDiff,acorr,lag,sig,a,b,t2,stat] = correlate(x1,x2,verbose,norm)
%Break out inputs
y1 = x1(:,2);
y2 = x2(:,2);

t1 = x1(:,3);
t2 = x2(:,3);

ub = min(t1(end),t2(end));

temp = (t2<=ub);
t2 = t2(temp);
y2 = y2(temp);
temp = (t1<=ub);
t1 = t1(temp);
y1 = y1(temp);

if isempty(t2)==1 || isempty(t1)==1
    error('times do not overlap')
    return
end



%Interpolate/Trim
y1_interp = interp1(t1,y1,t2);

if verbose
    figure
    plot(t2,y1_interp)
    hold on
    plot(t2,y2)
end

y1_interp(isnan(y1_interp))=[];
y2 = y2(1:length(y1_interp));
t2 = t2(1:length(y1_interp));

stat(1,1) = mean(y2);
stat(1,2) = std(y2);
stat(2,1) = mean(y1_interp);
stat(2,2)= std(y1_interp);

if norm
    %Normalize
    b=(y2-mean(y2))./std(y2);
    a=(y1_interp-mean(y1_interp))./std(y1_interp);
else
    
    b = y2;%;-mean(y2);
    a = y1_interp;%_interp-mean(y1_interp);
    
end

%XCorr
[acorr,lag] = xcorr(a,b);

%Results
[~,I] = max(abs(acorr));
lagDiff = lag(I);

n = length(y2);
k = lagDiff;
sig = 2/sqrt(n-k);
end