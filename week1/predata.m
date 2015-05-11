%%
T = readtable('sp500.csv');
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Index;
myfts = fints(dates, data);

%%
periods = {'12/30/1927::11/08/1929', '01/01/2007::12/31/2008','05/01/2013::04/30/2015','01/01/2010::12/31/2011'};
fts_window  = myfts(periods{1});

%%
y = fts2mat(fts_window.series1);
%prev = y(1);
%y(1) = [];
%logret = [log(y(1)) - log(prev); diff(log(y))];

LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
index = y;

%%
trendType = {'poly1', 'exp1'};

f=fit(x,y, trendType{2});
figure
plot(f, x, y);
y = y - f(x);

%%
f = fit(x,y,'sin2','Upper',[Inf 1 Inf Inf 1 Inf],'Robust','on');
figure
plot(f, x, y);
y = y - f(x);


