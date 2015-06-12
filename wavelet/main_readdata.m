defaultDates = '01/01/2007::01/20/2015';
eemusDates = '03/16/2011::01/20/2015';
fxDates = '01/07/2008::01/20/2015';
indexNames = {'spx','ukx','nky','dax','eemus','jpyusd','gbpusd','eurusd'};
imvNames = {'vix.mat','vftse.mat','vnky.mat','v1x.mat','vxeem.mat','jyvix.mat','bpvix.mat','euvix.mat'};

dataset = char(indexNames(dataindex));
imvfile = char(imvNames(dataindex));
dataPeriod = defaultDates;
if dataindex == 5
    dataPeriod = eemusDates;
elseif dataindex >= 6
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