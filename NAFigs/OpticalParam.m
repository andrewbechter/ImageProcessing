% Fiber optical parameters
%% Custom color lists, yo
d = get(groot,'DefaultAxesColorOrder');
for ii = 1:7
    colors{ii}=d(ii,:);
end
colors{8}= [0.175 0.175 0.175];
colors{9}= colors{2};
colors{10}= colors{3};

% % refractive index check
% lambda = [0.980,1.050];
% temp1=(0.6961663*lambda.^2)./(lambda.^2-0.0684043^2)+(0.4079426*lambda.^2)./(lambda.^2-0.1162414^2)+(0.8974794*lambda.^2)./(lambda.^2-9.896161^2);
% n = sqrt(temp1+1);
% figure;plot(n)

% ncl = 1.4507;
% nc = normrnd(1.45745,mean([9.0000e-04,9.99e-04]),1,ntracks);
% NA=sqrt(nc.^2-ncl.^2); %Numerical aperture of the fiber 1.4565-1.4584;


ntracks=500;
NA = 0.14;
nc = normrnd(0.14,0.003,1,ntracks);
% a= normrnd(4.75/2,0.35/2,1,ntracks); % physical core radius of the fiber in microns 4.4-5.1 microns
a= normrnd(4.82/2,0.35/2,1,ntracks); % physical core radius of the fiber in microns 4.4-5.1 microns
lambda = [0.950,0.980,1.05];

figure
for ii = 1:length(lambda)
    
    z = (0.1:1000:70000);
    V=(2*pi*a.*NA)./lambda(ii); %V number
    w0 =a.*(0.65+1.619./(V.^(1.5))+2.879./(V.^6)-(0.016+1.561./(V.^7))); % petermann_II MFR formula
    zr = pi*(w0.^2)./lambda(ii);
    temp = sqrt(1+((z'*(1./zr)).^2));
    w0 = ones(length(z),1)*w0;
    w(:,:,ii) = temp .*w0;
    mw(:,ii) = mean(w(:,:,ii)');
    sw(:,ii) = std(w(:,:,ii)');
    psigw(:,ii) = mw(:,ii)+sw(:,ii);
    nsigw(:,ii) = mw(:,ii)-sw(:,ii);
    plot(z,w(:,:,ii),'color',[colors{ii} 0.05],'linewidth',2)
    hold on
    
    %fitting
    lin= 'a1*x';
    [f,gof]=fit(z',mw(:,ii),lin,'Start',[0]);
    p = coeffvalues(f)
    
    [f,gof]=fit(z',psigw(:,ii),lin,'Start',[0]);
    ppsig = coeffvalues(f);
    
    [f,gof]=fit(z',nsigw(:,ii),lin,'Start',[0]);
    pnsig = coeffvalues(f);
    
    
%     y = f(z');
%     plot(z',mw(:,ii)','.','Color',colors{ii})
%     h(ii)=plot(z,y,'Color',colors{ii});
xlabel('Distance (\mum)')
ylabel('PSF radius (2\sigma) in \mum')
box on
grid on
% l=legend(h,{'950nm','1050nm','950nm 1','980nm 1','1050nm 1','950nm 2','980nm 2','1050nm 2','980nm 3','1050nm 3'});%,'1064nm'});
l=legend(h,{'950nm','980nm','1050nm'});%,'1064nm'});

ax = gca;
ax.GridAlpha = 0.1;
ax.GridLineStyle ='--';
ax.LineWidth = 0.75;
l.Color = 'none';
l.Location = 'Best';
    
end
% 
figure
hold on
% plot(z,psigw(:,1))
% plot(z,nsigw(:,1))
plot(z,mw)

xlabel('Distance (\mum)')
ylabel('PSF radius (2\sigma) in \mum')
box on
grid on
% l=legend(h,{'950nm','1050nm','950nm 1','980nm 1','1050nm 1','950nm 2','980nm 2','1050nm 2','980nm 3','1050nm 3'});%,'1064nm'});
l=legend(h,{'950nm','980nm','1050nm'});%,'1064nm'});

ax = gca;
ax.GridAlpha = 0.1;
ax.GridLineStyle ='--';
ax.LineWidth = 0.75;
l.Color = 'none';
l.Location = 'Best';


