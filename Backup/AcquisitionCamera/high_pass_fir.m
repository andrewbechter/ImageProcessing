function [high_data] = high_pass_fir(p,data)
Fstop = p(1,1);
Fpass = p(1,2);
Astop = p(1,3);
Apass = p(1,4);
Fs = p(1,5);

d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');
fvtool(d)
high_data = filter(d,data);
end
