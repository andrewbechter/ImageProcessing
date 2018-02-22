%omega progagation for theretical NA of SMF
clear
%Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end
colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};
clear d
nc = 1.4568;%1.4565-1.4584; % This value gives the right "NA" from the spec sheet 5.8
ncl = 1.4507; 
figure()
hold on
ii = 1;

for lambda = [0.882,0.950,0.980,1.050,1.31]
na=sqrt(nc^2-ncl^2); %Numerical aperture of the fiber
na = 0.136;
a=4.95/2; % physical core radius of the fiber in microns
V=(2*pi*a*na)/lambda; %V number
w0 =a.*(0.65+1.619./(V.^(1.5))+2.879./(V.^6)-(0.016+1.561./(V.^7))); % petermann_II MFR formula
waist(ii) = w0;

z = (30000:5:70000); %propagation distnace
zr = pi*(w0^2)/lambda;
w = w0*sqrt(1+(z./zr).^2);
% R = x.*(1+(zr./z).^2); % Radius of curvature (wavefront) 
w1 = w-w(1);
z1 = z-z(1);

w = w1;
z = z1;

hold on
plot(z,w,'.','Color',colors{ii})
lin= 'a1*x';
[f,gof]=fit(z',w',lin,'Start',[0]);
y = f(z);
h(ii)=plot(z',y','Color',colors{ii});
p(ii,1) = coeffvalues(f);
disp([lambda,p(1,1),w0,V])
ii = ii+1;
end

xlabel('Distance (\mum)')
ylabel('PSF radius (2\sigma) in \mum')
box on
grid on
l=legend(h,{'882nm','950nm','980nm','1050nm','1310nm'});%,'1064nm'});
ax = gca;
ax.GridAlpha = 0.1;
ax.GridLineStyle ='--';
ax.LineWidth = 0.75;
l.Color = 'none';
l.Location = 'Best';

% Create textbox
annotation('textbox',[0.2385 0.6 0.35 0.25],...
    'String',{['NA =',num2str(p(1,1)),' @\lambda = 882nm'],['NA =',num2str(p(1,1)),' @\lambda = 950nm'],['NA =', num2str(p(2,1)),' @\lambda = 980nm'],['NA =', num2str(p(3,1)),' @\lambda = 1050nm'],['NA =',num2str(p(5,1)),' @\lambda = 1310nm']});








