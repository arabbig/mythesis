
T = readtable('sp500.csv');
dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
data = T.Index;
myfts = fints(dates, data);

%%
str1 = '12/30/1927::12/31/1929';
str2 = '01/01/2008::12/31/2008';
fts_window  = myfts(str2);


%%
y = fts2mat(fts_window.series1);
prev = y(1);
y(1) = [];
logret = [log(y(1)) - log(prev); diff(log(y))];

LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
index = y;

%%
f=fit(x,y,'poly1');

%%
figure
plot(f, x, y);
y = y - f(x);

%%
f=fit(x,y,'exp1');

%%
figure
plot(f, x, y);
y = y - f(x);

%%
f = fit(x, y, 'poly2','Robust','on');
%%
figure;
plot(f, x, y);
y = y - f(x);
%%
f=fit(x, y, 'sin1','Robust','on');
%%
figure;
plot(f, x, y);
y = y - f(x);

%% 
%See https://www.kevinsheppard.com/images/9/95/MFE_Toolbox_Documentation.pdf

[trend, cyclic] = bkfilter(y,6,32);
y = cyclic;

%%
figure
subplot(3,1,1)
coefs = cwt(y,1:128,'meyr','plot');
title('Continuous Transform, absolute coefficients');
set(gca,'yticklabel',[]);
ylabel('Scale');

subplot(3,1,2)
plot(x,index,'b')
title('Index')

%%
coefs_by_day = num2cell(coefs,1);
coefs_med = cellfun(@median, coefs_by_day);
coefs_mean = cellfun(@mean, coefs_by_day);
qfunc = @(q) @(x) quantile(x,q);
ninethy = qfunc(.90);
coefs_tail = cellfun(ninethy, coefs_by_day);
coefs_max = cellfun(@max, coefs_by_day);
subplot(3,1,3);
plot(x,coefs_mean,'r');
hold on
plot(x,coefs_med,'g');
hold on
plot(x,coefs_tail,'b');
hold on
plot(x,coefs_max,'c');

