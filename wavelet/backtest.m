
    
%% Read whole period data
T = readtable([dataset,'.csv']);
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Value;
fts = fints(dates, data);

%% Calculating CWT, Zetas
fts_window  = fts(backtestPeriod);
N = length(fts_window);
windowSize = 252;
zetas_lookback1 = ones(N-windowSize, 1);
coefs_lookback1 = ones(32,N-windowSize);
zetas_lookback15 = ones(N-windowSize, 1);
coefs_lookback15 = ones(32,N-windowSize);

RANGE = (windowSize + 1) : N ;
for date = RANGE
    coefs = calcwt(fts_window((date - windowSize):(date - 1)));
    cf_lastday = coefs(:,end);
    cf_last15day = coefs(:,end-15);
    zetas_lookback1(date-windowSize) = calzeta(cf_lastday,coefs,'global_mean');
    coefs_lookback1(:,date-windowSize) = cf_lastday;
    zetas_lookback15(date-windowSize) = calzeta(cf_last15day,coefs,'global_mean');
    coefs_lookback15(:,date-windowSize) = cf_last15day;
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

axeAct = subplot(3,1,1);
image(wcodemat(abs(coefs_lookahead),240,'mat',1))
set(axeAct,'YDir','normal')


axeAct = subplot(3,1,2);
image(wcodemat(abs(coefs_lookback1),240,'mat',1))
set(axeAct,'YDir','normal')
%%
coefs_cell = num2cell(coefs_lookahead,1);
calzt = @(c) calzeta(c,coefs_lookahead,'global_mean');
zetas_lookahead = cellfun(calzt,coefs_cell);

subplot(3,1,1);
% Set axes http://uk.mathworks.com/matlabcentral/newsreader/view_thread/148694
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_lookahead);
datetick(hxs2(1));
datetick(hxs2(2));

subplot(3,1,2);
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_lookback1);
datetick(hxs2(1));
datetick(hxs2(2));

subplot(3,1,3);
[hxs2,~,~]=plotyy(testTS.dates,fts2mat(testTS.series1),testTS.dates,zetas_lookback15);
datetick(hxs2(1));
datetick(hxs2(2));

%%
figure;
%subplot(3,1,2)
caldrawdown(testTS,zetas_lookahead,0.0);
%subplot(3,1,2);
%caldrawdown(testTS,zetas2,0.0);

%%
RANGE2 = 20:length(testTS);
logret = [0; diff(log(fts2mat(testTS.series1)))];
monthvols = sqrt(conv(logret.^2,ones(20,1))*252/20);

subplot(2,1,1);
[hxs,~,~]=plotyy(testTS.dates(RANGE2),zetas_lookback1(RANGE2),testTS.dates(RANGE2),monthvols(RANGE2));
datetick(hxs(1));
datetick(hxs(2));

subplot(2,1,2);
[hxs1,~,~]=plotyy(testTS.dates(RANGE2),fts2mat(testTS.series1(RANGE2)),testTS.dates(RANGE2),zetas_lookback1(RANGE2));
datetick(hxs1(1));
datetick(hxs1(2));
