trendType = {'poly1', 'exp1'};
%%
T = readtable('sp500.csv');
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Index;
myfts = fints(dates, data);
%%
periods = {'12/30/1927::11/08/1929', '01/01/2007::12/31/2008','05/01/2013::04/30/2015','01/01/2010::12/31/2011'};
fts_window  = myfts(periods{1});

%%
y = fts2mat(fts_window.series1);
prev = y(1);
y(1) = [];
logret = [log(y(1)) - log(prev); diff(log(y))];

LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
index = y;

%%
f=fit(x,y, trendType{2});
figure
plot(f, x, y);
y = y - f(x);

%%
f = fit(x,y,'sin2','Upper',[Inf 1 Inf Inf 1 Inf],'Robust','on');
figure
plot(f, x, y);
y = y - f(x);

%% 
% See https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf

% [trend, cyclic] = bkfilter(y,6,32);
% y = cyclic;

%%
figure
subplot(3,1,1)
coefs = cwt(y,1:32,'meyr','plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');
subplot(3,1,2)
plot(fts_window);
title('Index')

%%

coefs_by_day = num2cell(coefs,1);
morethan = @(tresh) @(x) sum(x>tresh)/length(x);
eta = morethan(0);
coefs_eta = cellfun(eta, coefs_by_day);
hold on
subplot(3,1,3);
plot(x,coefs_eta)

%%

pfunc = @(q) @(x) prctile(x,q);
ninethyP = pfunc(90);
tenthP = pfunc(10);

tenthP_ccfs = cellfun(tenthP, coefs_by_day);
median_ccfs = cellfun(@median, coefs_by_day);
ninethyP_ccfs = cellfun(ninethyP, coefs_by_day);

subplot(3,1,3);
plot(x, median_ccfs,'r');
hold on
plot(x,tenthP_ccfs,'g');
hold on
plot(x,ninethyP_ccfs,'b');
