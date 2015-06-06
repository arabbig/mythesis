defaultDates = '01/01/2007::01/20/2014';
eemusDates = '03/16/2011::01/20/2014';
fxDates = '01/07/2008::01/20/2014';
indexNames = {'spx','eemus','ukx','nky','dax','jpyusd','gbpusd','eurusd'};
imvNames = {'vix.mat','vxeem.mat','vftse.mat','vnky.mat','v1x.mat','jyvix.mat','bpvix.mat','euvix.mat'};

dataset = char(indexNames(i));
imvfile = char(imvNames(i));
dataPeriod = defaultDates;
if i == 2
    dataPeriod = eemusDates;
elseif i >= 6
    dataPeriod = fxDates;
end
%%
load(imvfile);
fprintf('data drawn from = %s \n', dataPeriod);
fprintf('Underlying = %s \n', dataset);
fprintf('Implied vol = %s \n', imvfile(1:end-4));

%% Get windowed time series 
T = readtable([dataset,'.csv']);
fts =  createFTS(T);
ivol_fts = createFTS(imvTab);

fts_window  = fts(dataPeriod);

%% Calculating CWT, Zetas
N = min(length(fts_window));
windowSize = 252;

zetas_delay5 = ones(N-windowSize, 1);
zetas_mean10 = ones(N-windowSize, 1);
%zetas_max5 = ones(N-windowSize, 1);
%zetas_max10 = ones(N-windowSize, 1);
RANGE = (windowSize + 1) : N ;
for date = RANGE 
    coefs = calcwt(fts_window((date - windowSize):(date-1)));
    tresh = calthreshold(coefs,'global_mean');
    zs = ones(10,1);
    for i = 1:10
        cf = coefs(:,end-i+1);
        zs(i) = calzeta2(cf, tresh); 
    end
    zetas_delay5(date-windowSize) = zs(4);
    zetas_mean10(date-windowSize) = mean(zs);    
end
%%
testTS = fts_window(RANGE);
logret = diff(log(fts2mat(testTS.series1)));
realVolMAT = sqrt(conv(logret.^2,ones(20,1))*252/20)*100;
% align it with their corresponding testTS.dates
% check
% realized(t) = sqrt(mean(diff(log(fts2mat(testTS.series1(t-20:t)))).^2)*252)*100
% == nxtmnthVols(t) except precision loss
realVolMAT(end-18:end) = [];
realVolMAT(1:19) = NaN;
realVolMAT = [NaN; realVolMAT]; 
realVolTS = fints(testTS.dates,realVolMAT);
ivol_window  = ivol_fts(dataPeriod);

testTS = chfield(testTS, 'series1', 'index');
realVolTS = chfield(realVolTS, 'series1', 'realVol');
ivol_window = chfield(ivol_window, 'series1', 'impVol');
allTS = merge(testTS,realVolTS);
allTS = merge(allTS, ivol_window); % merge to make ivol aligned with index, in case there are missing values
impVolTS_lag20 = lagts(extfield(allTS,'impVol'),20);
spreadTS = chfield(extfield(allTS,'realVol'),'realVol','spread') - chfield(lagts(extfield(allTS,'impVol'),20),'impVol','spread');
spreadTS = spreadTS(datestr(testTS.dates));

%missing data in BPVIXIndex after 29 May 2009

%%
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, zetas_mean10);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('Mean 10 days')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(testTS.dates, sign(fts2mat(spreadTS.spread)), testTS.dates, zetas_mean10);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
pair = [zetas_mean10, fts2mat(spreadTS.spread)];
pair(any(isnan(pair),2),:) = [];
crosscorr(pair(:,1),pair(:,2),100);

%%
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(testTS.dates, fts2mat(spreadTS.spread), testTS.dates, zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('Delay 5 days')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(testTS.dates, sign(fts2mat(spreadTS.spread)), testTS.dates, zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
pair = [zetas_delay5, fts2mat(spreadTS.spread)];
pair(any(isnan(pair),2),:) = [];
crosscorr(pair(:,1),pair(:,2),100);