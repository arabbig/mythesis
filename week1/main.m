asset(1).name = 'spx';
asset(1).periods = {'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};

asset(2).name = 'dax';
asset(2).periods = {'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};

asset(3).name = 'nky';
asset(3).periods = {'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};
%http://en.wikipedia.org/wiki/Japanese_asset_price_bubble

asset(4).name = 'ukx';
asset(4).periods = {'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};

asset(5).name = 'eurusd';
asset(5).periods = {'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};

asset(6).name = 'gbpusd';
asset(6).periods = {'01/01/2008::12/31/2009','01/01/2009::12/31/2010'};

asset(7).name = 'jpyusd';
asset(7).periods = {'01/01/2012::12/31/2013','01/01/2013::12/31/2014'};

%%
assetID = 6;
ast = asset(assetID);
dataset=ast.name;
periods= ast.periods;
periodID = 1;

%% Read whole period data
T = readtable([dataset,'.csv']);
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Value;
fts = fints(dates, data);

%% Calculating CWT, Zetas
thisPeriod = periods{periodID};
disp(thisPeriod);
fts_window  = fts(thisPeriod);
coefs_zeta = calcwt(fts_window);

%% Calculating drawdown
minzet = 0.4;
fid = fopen(['output/',dataset,'_',thisPeriod(7:10),'_',thisPeriod(19:22),'_min',num2str(minzet*10),'.txt'],'wt');
caldrawdown(fts_window, coefs_zeta, minzet, fid);
fclose(fid);
