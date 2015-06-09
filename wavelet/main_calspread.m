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

testTS = chfield(testTS, 'series1', 'index');
realVolTS = chfield(realVolTS, 'series1', 'realVol');
ivol_window = chfield(ivol_window, 'series1', 'impVol');
allTS = merge(testTS,realVolTS);
allTS = merge(allTS, ivol_window); % merge to make ivol aligned with index, in case there are missing values
impVolTS_lag20 = lagts(extfield(allTS,'impVol'),20);
spreadTS = chfield(extfield(allTS,'realVol'),'realVol','spread') - chfield(lagts(extfield(allTS,'impVol'),20),'impVol','spread');
spreadTS = spreadTS(datestr(testTS.dates));

%missing data in BPVIXIndex after 29 May 2009