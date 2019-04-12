%   Thanks to R. Pawlowicz (IOS) rich@ios.bc.ca for 'contours.m' and 
%   'clabel.m/inline_labels' so that contour now works with parametric
%   surfaces and inline contour labels.
t=0:0.1:7;
y=sin(t);
z=cos(t);
figure
plot(t,t,'r',t,y,'g',t,z,'b',t,y+z,'k');
Strings={'line','sin','cos','sin+cos'};
LabelLine(gca,Strings);

function LabelLine(h_axes,LabelStrings)
h_line=findobj(h_axes,'type','line');
n=length(h_line);
for k=1:n;
    XData=get(h_line(k),'XData');
    YData=get(h_line(k),'YData');
    Ind=round(length(XData)/2);
    text(XData(Ind),YData(Ind),LabelStrings{n-k+1});
end
end