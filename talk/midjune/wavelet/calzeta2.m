function zeta = calzeta2(coef,tresh)
morethan = @(tresh) @(x) sum(x>tresh)/length(x);
zetafun = morethan(tresh);
zeta = zetafun(coef);
end