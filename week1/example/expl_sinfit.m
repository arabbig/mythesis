% Examples are from
% http://uk.mathworks.com/help/curvefit/least-squares-fitting.html

xdata = (.25*pi:0.1:6*pi)';
y0 = sin(xdata);

% Response-dependent Gaussian noise
gnoise = y0.*randn(size(y0));

% Salt-and-pepper noise
spnoise = zeros(size(y0));
p = randperm(length(y0));
sppoints = p(1:round(length(p)/5));
spnoise(sppoints) = 5*sign(y0(sppoints));

ydata = y0 + gnoise + spnoise;

f = fittype('a*sin(b*x)');

[fit1,gof,fitinfo] = fit(xdata,ydata,f,'StartPoint',[1 1]);
residuals = fitinfo.residuals;
I = abs( residuals) > 1.5 * std( residuals );
outliers = excludedata(xdata,ydata,'indices',I);
fit2 = fit(xdata,ydata,f,'StartPoint',[1 1],...
           'Exclude',outliers);
fit3 = fit(xdata,ydata,f,'StartPoint',[1 1],'Robust','on');
hold on
plot(xdata,y0,'g.-')
plot(fit1,'r-',xdata,ydata,'k.',outliers,'m*')
plot(fit2,'c--')
plot(fit3,'b:')