
    
%% Read whole period data
T = readtable([dataset,'.csv']);
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Value;
fts = fints(dates, data);

%% Calculating CWT, Zetas
fts_window  = fts(backtestPeriod);
N = length(fts_window);
windowSize = 252;
zetas1 = ones(N-windowSize, 1);
%zetas2 = ones(N-windowSize, 1);
%ccfs1 = ones(32,N-windowSize);
%ccfs2 = ones(32,N-windowSize);

RANGE = (windowSize + 1) : N ;
for date = RANGE
    [zeta1 ccf1] = calcwt(fts_window((date - windowSize):(date - 1)),'global_mean');
    %[zeta2 ccf2] = calcwt(fts_window((date - windowSize):(date - 1)),'local_mean');
    zetas1(date-windowSize) = zeta1;
    %zetas2(date-windowSize) = zeta2;
    %ccfs1(:,date-windowSize) = ccf1;
    %ccfs2(:,date-windowSize) = ccf2;
    %disp(date-windowSize);
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
figure;
subplot(3,1,1);
coefs = cwt(y,1:32,'meyr','plot');

%%
subplot(3,1,2)
caldrawdown(testTS,zetas1,0.0);
%subplot(3,1,2);
%caldrawdown(testTS,zetas2,0.0);

%%
RANGE2 = 20:length(testTS);
logret = [0; diff(log(fts2mat(testTS.series1)))];
monthvols = sqrt(conv(logret.^2,ones(20,1))*252/20);

subplot(3,1,2);
[hxs,~,~]=plotyy(testTS.dates(RANGE2),zetas1(RANGE2),testTS.dates(RANGE2),monthvols(RANGE2));
datetick(hxs(1));
datetick(hxs(2));

subplot(3,1,3);
[hxs1,~,~]=plotyy(testTS.dates(RANGE2),fts2mat(testTS.series1(RANGE2)),testTS.dates(RANGE2),monthvols(RANGE2));
datetick(hxs1(1));
datetick(hxs1(2));
