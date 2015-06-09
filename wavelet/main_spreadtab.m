zetaD5TS = lagts(fints(testTS.dates,zetas_delay5,'zeta'),20);
zetaM5TS = lagts(fints(testTS.dates,zetas_mean10,'zeta'),20);

fprintf('Days during testing = %d\n',length(zetaD5TS));
fprintf('   X      days     AvgSpread    days    AvgSpread \n');
fprintf('        (delay5)    (delay5)  (mean10)   (mean10) \n');
for X = linspace(0,0.9,10)
    ii  = fts2mat(zetaD5TS.zeta) > X;
    jj  = fts2mat(zetaM5TS.zeta) > X;
    fprintf('%6.3f%9.d%11.3f%10.1d%11.3f\n', X, sum(ii),nanmean(-spreadTS(ii)),sum(jj)+0, nanmean(-spreadTS(jj)));
end
