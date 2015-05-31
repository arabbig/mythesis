spec2 = ugarchspec(mean.model=list(armaOrder=c(2,2)
fit  = ugarchfit(data = y, spec= spec2)
grchvol = as.numeric(sigma(fit));
plot(y,ty='l')
lines(grchvol,col='red')
lines(sqrt(mvol),col='blue')
lines(exp(x/2),col='green')

