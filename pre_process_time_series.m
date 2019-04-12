function [x1,x2] = pre_process_time_series(x1,x2,Fs)

offset = (x2(1,1)-x1(1,1))*24*3600;  

if abs(offset) > 1000
    x2(:,1) = x2(:,1)+7/24;
    offset = (x2(1,1)-x1(1,1))*24*3600;  
end

step1 = 10e-3;
step2 = 1./Fs;

t1 = (0:1:length(x1(:,1))-1)*step1;
t2 = (0:1:length(x2(:,1))-1)*step2;

if offset > 1
    t2 = t2+offset;
else
    t1 = t1+abs(offset);
end

x1(:,3) = t1;
x2(:,3) = t2;
end