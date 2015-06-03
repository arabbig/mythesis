dataset = 'spx';
backtestPeriod = '01/01/2007::01/20/2011';
    
%% Read whole period data
T = readtable([dataset,'.csv']);
fts =  createFTS(T);

%% Calculating CWT, Zetas
fts_window  = fts(backtestPeriod);
N = length(fts_window) - 15;
windowSize = 252;

zetas_pad0 = ones(N-windowSize, 1);
coefs_pad0 = ones(32,N-windowSize);
zetas_pad5 = ones(N-windowSize, 1);
coefs_pad5 = ones(32,N-windowSize);
zetas_pad16 = ones(N-windowSize, 1);
coefs_pad16 = ones(32,N-windowSize);

RANGE = (windowSize + 1) : N ;
for date = RANGE
    coefs = calcwt(fts_window((date - windowSize):(date)));
    cf_today = coefs(:,end);
    zetas_pad0(date-windowSize) = calzeta(cf_today,coefs,'global_mean');
    coefs_pad0(:,date-windowSize) = cf_today;

    coefs = calcwt(fts_window((date - windowSize):(date+5)));
    cf_today_pad5 = coefs(:,end-5);
    zetas_pad5(date-windowSize) = calzeta(cf_today_pad5,coefs,'global_mean');
    coefs_pad5(:,date-windowSize) = cf_today_pad5;    
    
    coefs = calcwt(fts_window((date - windowSize):(date+15)));
    cf_today_pad16 = coefs(:,end-15);
    zetas_pad16(date-windowSize) = calzeta(cf_today_pad16,coefs,'global_mean');
    coefs_pad16(:,date-windowSize) = cf_today_pad16;
end
%%
testTS = fts_window(RANGE);

%%
y = fts2mat(testTS);
LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
trendType = {'poly1', 'exp1'};
f=fit(x,y, trendType{1});
y = y - f(x);
coefs_lookahead = cwt(y,1:32,'meyr');

%%
figure;
colormap(pink(240))

axeAct = subplot(4,1,1);
image(wcodemat(abs(coefs_lookahead),240,'mat',1))
set(axeAct,'YDir','normal')
title('Whole period CWT')

axeAct = subplot(4,1,2);
image(wcodemat(abs(coefs_pad16),240,'mat',1))
set(axeAct,'YDir','normal')
title('Backtest, 16 days padding')

axeAct = subplot(4,1,3);
image(wcodemat(abs(coefs_pad5),240,'mat',1))
set(axeAct,'YDir','normal')
title('Backtest, 5 days padding')

axeAct = subplot(4,1,4);
image(wcodemat(abs(coefs_pad0),240,'mat',1))
set(axeAct,'YDir','normal')
title('Backtest, no padding')
%%
figure;
coefs_cell = num2cell(coefs_lookahead,1);
calzt = @(c) calzeta(c,coefs_lookahead,'global_mean');
zetas_lookahead = cellfun(calzt,coefs_cell);

subplot(4,1,1);
% Set axes http://uk.mathworks.com/matlabcentral/newsreader/view_thread/148694
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_lookahead);
datetick(hxs2(1));
datetick(hxs2(2));
title('Whole period Zeta')

subplot(4,1,2);
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_pad16);
datetick(hxs2(1));
datetick(hxs2(2));
title('Backtest, 16 days padding')

subplot(4,1,3);
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_pad5);
datetick(hxs2(1));
datetick(hxs2(2));
title('Backtest, 5 days padding')

subplot(4,1,4);
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_pad0);
datetick(hxs2(1));
datetick(hxs2(2));
title('Backtest, No padding')

%%
%figure;
%subplot(3,1,2)
%caldrawdown(testTS,zetas_lookahead,0.0);
%subplot(3,1,2);
%caldrawdown(testTS,zetas2,0.0);

