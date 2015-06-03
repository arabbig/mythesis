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
ivol_window  = ivol_fts(dataPeriod);

%% Calculating CWT, Zetas
N = min(length(fts_window));
windowSize = 252;

zetas_delay0 = ones(N-windowSize, 1);
zetas_delay5 = ones(N-windowSize, 1);

RANGE = (windowSize + 1) : N ;
for date = RANGE 
    coefs = calcwt(fts_window((date - windowSize):(date)));
    cf_today = coefs(:,end);
    zetas_delay0(date-windowSize) = calzeta(cf_today, coefs,'global_mean');   
    cf_past5day = coefs(:,end-5);
    zetas_delay5(date-windowSize) = calzeta(cf_past5day, coefs,'global_mean');    
end
%%
testTS = fts_window(RANGE);
logret = diff(log(fts2mat(testTS.series1)));
nxtmnthVols = sqrt(conv(logret.^2,ones(20,1))*252/20)*100;
nxtmnthVols(end-18:end) = [];
nxtmnthVols(1:19) = [];
L = length(nxtmnthVols);
realVolTS = fints(testTS.dates(1:L)+20,nxtmnthVols);
[d,ia,ib]=intersect(realVolTS.dates, ivol_window.dates);
realVolTS = realVolTS(ia);
zetas_delay0 = zetas_delay0(ia);
zetas_delay5 = zetas_delay5(ia);
ivolTS = ivol_window(ib);

%missing data in BPVIXIndex after 29 May 2009

%%
loss = fts2mat(realVolTS) - fts2mat(ivolTS.series1);
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(realVolTS.dates,loss,realVolTS.dates,zetas_delay0);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('No delay')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(realVolTS.dates,sign(loss),realVolTS.dates,zetas_delay0);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
crosscorr(zetas_delay0,loss);
%%
figure;
subplot(3,1,1);
[hxs,h1,h2]=plotyy(realVolTS.dates,loss,realVolTS.dates,zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')
title('5 days delay')

subplot(3,1,2);
[hxs,h1,h2]=plotyy(realVolTS.dates,sign(loss),realVolTS.dates,zetas_delay5);
datetick(hxs(1));
datetick(hxs(2));
set(hxs(1),'YLim',[-2,2])
set(hxs(2),'YLim',[0,3])
legend([h2;h1],'zeta','spread')

subplot(3,1,3);
crosscorr(zetas_delay5,loss);