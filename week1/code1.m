% read data in
T = readtable('sp500.csv');

Index = table2array(T(:,{'Index'}));

figure
%subplot(2,1,2)
prc = Index(end-10*250:end)
L = length(prc)
TimeInYear = (0:(L-1))/250 + 2005;
plot(TimeInYear, prc)
title('Last 10 years')

%% 

yr = 3; %3 years back data look linear and cyclical
Y = Index(end-yr*250:end);
L = length(Y);

X = ones(L, 2);
X(:,2) = 0:(L-1);
%time normalized into yr
X(:,2) = X(:,2)/250;
beta = X\Y;


TimeInYear = X(:,2)
figure
subplot(2,1,1)
plot(TimeInYear,Y,'b')
title('Last 3 Years')
hold on
plot(TimeInYear,X*beta,'r')
subplot(2,1,2)
Y_det = Y - X*beta; % detrend
plot(TimeInYear, Y_det,'.b')
title('Last 3 Years - Detrend')

%% 
%See http://uk.mathworks.com/help/releases/R2014b/curvefit/sum-of-sine.html

%str = 'sin1'
str = 'sin6'

f=fit(TimeInYear,Y_det, str)
figure
subplot(2,1,1)
plot(f,TimeInYear,Y_det)
title(str)
subplot(2,1,2)
plot(f,TimeInYear,Y_det,'residuals')
title('Residuals')

%% 
%See https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf

[trend, cyclic] = bkfilter(Y_det,6,32)
figure
subplot(2,1,1)
plot(TimeInYear,trend,'r',TimeInYear,Y_det,'.b')
title('Baxter-King Filter')
subplot(2,1,2)
plot(TimeInYear,cyclic,'.b')
hline = refline([0 0])
hline.Color = 'r'
title('Residuals')

%% 
%See http://uk.mathworks.com/help/signal/ref/periodogram.html

[pxx,w] = periodogram(Y_det);
[pxx1,w1] = periodogram(cyclic);

figure
subplot(2,1,1)
plot(w/2/pi,pxx,'b')
hold on
plot(w1/2/pi,pxx1,'r')
xlabel('\omega / \pi')
title('Periodogram')
legend('Orginal','High-pass')


subplot(2,1,2)
plot(w/2/pi,pxx,'b')
hold on
plot(w1/2/pi,pxx1,'r')
xlabel('\omega / \pi')
legend('Orginal','High-pass')
axis([0 inf -2000,2000])

%%
figure
subplot(2,1,1)
ccfs=cwt(cyclic,1:128,'sym2','plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');
subplot(2,1,2)
plot(TimeInYear,Y,'b')
title('Original')

%%

