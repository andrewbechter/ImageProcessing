%contour test
load('/Volumes/Software/Simulator/TiptiltStudy/tt_v_fc.mat')
x = data(:,1);
y = x;
z = data(:,2);

F = scatteredInterpolant(x,y,z);
Z = F(x,y); %evaluate surface

imagesc(x,y,Z)

fit(data(:,1),data(:,2),'gauss1')



return


v = [15,10,5];

figure
[C,h] = contour(X,Y,Z,v,'color','k','linestyle','--','linewidth',1.5);
clabel(C,h)

