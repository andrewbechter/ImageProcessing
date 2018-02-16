function [smooth_data] = Moving_Average(windowSize,data)
%Moving Average Filter%
% windowSize = 1000;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smooth_data = filter(b,a,data);
end