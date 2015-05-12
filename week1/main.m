dataset='sp500';
periodNum = 2;

T = readtable([dataset,'.csv']);
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Index;
fts = fints(dates, data);
periods={'01/01/2007::12/31/2008','01/01/2010::12/31/2011'};
thisPeriod = periods{periodNum};
fts_window  = fts(thisPeriod);

coefs_zeta = calcwt(fts_window);

%%
fid = fopen(['output/',dataset,'_',thisPeriod(7:10),'_',thisPeriod(19:22),'.txt'],'wt');
caldrawdown(fts_window, coefs_zeta,0.2, fid);
fclose(fid);
