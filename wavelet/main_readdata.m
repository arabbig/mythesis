defaultDates = '01/01/2007::01/20/2015';
eemusDates = '03/16/2011::01/20/2015';
fxDates = '01/07/2008::01/20/2015';
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

%%
ivol_window  = ivol_fts(dataPeriod);
fts_window  = fts(dataPeriod);