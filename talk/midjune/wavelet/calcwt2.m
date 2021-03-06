function coefs=calcwt2(fts_window,padSize)
y = fts2mat(fts_window.series1);
LEN = length(y);
x = (0:LEN-1)'/(LEN-1);
trendType = {'poly1', 'exp1'};
f=fit(x,y, trendType{1});
y = y - f(x);
y(end+1:end+padSize) = y(end);
coefs = cwt(y,1:32,'meyr');
end