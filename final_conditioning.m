function [y1_trim,y2_trim,x1_trim,gamma] = final_conditioning(indDiff,acorr,lag,y1,y2,x,Fs,verbose)

timeDiff = indDiff/Fs;

if verbose
figure
plot(lag,acorr)
end

x1 = x+timeDiff;

lb = max(x1(1),x(1));
ub = min(x1(end),x(end));

x1_trim = (lb:1/Fs:ub);

y1_trim = interp1(x,y1,x1_trim);
y2_trim = interp1(x1,y2,x1_trim);

gamma = corrcoef(y2_trim,y1_trim);
end