function zeta = calzeta(coef,coefs,str)
thr = calthreshold(coefs,str);
morethan = @(tresh) @(x) sum(x>tresh)/length(x);
zetafun = morethan(thr);
zeta = zetafun(coef);
end